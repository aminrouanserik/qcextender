import numpy as np
from typing import Iterable, Sequence, Self
from scipy.interpolate import make_interp_spline
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.types import timeseries as ts
from pycbc.filter.matchedfilter import match as cbcmatch
from qcextender.metadata import Metadata
from qcextender.basewaveform import BaseWaveform
from qcextender.models import lal_mode
from qcextender.functions import spherical_harmonics


class Waveform(BaseWaveform):
    """Waveform class with which all calculations can be done, stores multi-modal (time domain) waveforms."""

    def __init__(
        self, strain: np.ndarray, time: np.ndarray, metadata: Metadata
    ) -> None:
        """Initializes class containing waveform data.

        Args:
            strain (np.ndarray): Stacked array of multi-modal wave strains.
            time (np.ndarray): Time array, should be the same length as component strain arrays.
            metadata (dict): Metadata belonging to the generated or requested waveform.
        """
        super().__init__(strain, time, metadata)

    @classmethod
    def from_model(
        cls, approximant: str, modes: Iterable[Sequence[int]] = [(2, 2)], **kwargs
    ) -> Self:
        """Generates a multi-modal time-domain waveform from a specified model using PyCBC.

        Args:
            approximant: PyCBC included waveform model.
            modes: Modes to include, each as a (l, m) pair.

        Returns:
            Waveform: Object with stacked modes as complex array.
        """
        total_mass = kwargs["mass1"] + kwargs["mass2"]

        q = kwargs["mass1"] / kwargs["mass2"]
        if q < 1:
            q = 1 / q

        kwargs.update(
            library="lalsimulation",
            q=q,
            approximant=approximant,
            modes=list(modes),
            total_mass=total_mass,
        )
        metadata = cls._kwargs_to_metadata(kwargs)

        single_mode_strain = []
        if approximant in ["IMRPhenomD", "SEOBNRv4"]:
            for mode in modes:
                time, strain = lal_mode(
                    approximant,
                    kwargs["mass1"],
                    kwargs["mass2"],
                    metadata["spin1"],
                    metadata["spin2"],
                    metadata["distance"],
                    metadata["coa_phase"],
                    metadata["delta_t"],
                    metadata["f_lower"],
                    kwargs["f_ref"],
                    mode,
                )
            single_mode_strain.append(strain)

        multi_mode_strain = np.vstack(single_mode_strain)
        time = cls._align(single_mode_strain[0], time)

        return cls(
            multi_mode_strain,
            time,
            metadata,
        )

    def match(
        self,
        waveform: Self,
        f_lower: float | None = None,
        psd: str = "aLIGOZeroDetHighPower",
    ) -> float:
        """Calculates the match between self and one other waveform. Note, only the real parts are used.

        Args:
            waveform (Self): A waveform to match with.
            f_lower (float | None, optional): Lower cut-off for the match. Defaults to None, which then takes the highest of the two waveforms from the metadata.
            psd (str, optional): The PSD to use in the match. Defaults to "aLIGOZeroDetHighPower".

        Returns:
            float: The match of the real part of two modes.
        """
        delta_t = max(self.metadata.delta_t, waveform.metadata.delta_t)
        wf1_time = np.arange(self.time[0], self.time[-1], delta_t)
        wf2_time = np.arange(waveform.time[0], waveform.time[-1], delta_t)

        wf1_strain = 0
        for mode in self.metadata.modes:
            single_mode = make_interp_spline(self.time, self[mode])(wf1_time)
            single_minus_mode = make_interp_spline(self.time, self[mode[0], -mode[1]])(
                wf1_time
            )
            wf1_strain += single_mode * spherical_harmonics(
                mode[0], mode[1], self.metadata.inclination, self.metadata.coa_phase
            ) + single_minus_mode * spherical_harmonics(
                mode[0], -mode[1], self.metadata.inclination, self.metadata.coa_phase
            )

        wf2_strain = 0
        for mode in waveform.metadata.modes:
            single_mode = make_interp_spline(waveform.time, waveform[mode])(wf2_time)
            single_minus_mode = make_interp_spline(
                waveform.time, waveform[mode[0], -mode[1]]
            )(wf2_time)
            wf2_strain += single_mode * spherical_harmonics(
                mode[0],
                mode[1],
                waveform.metadata.inclination,
                waveform.metadata.coa_phase,
            ) + single_minus_mode * spherical_harmonics(
                mode[0],
                -mode[1],
                waveform.metadata.inclination,
                waveform.metadata.coa_phase,
            )

        wf1 = ts.TimeSeries(wf1_strain.real, delta_t=delta_t)
        wf2 = ts.TimeSeries(wf2_strain.real, delta_t=delta_t)

        if f_lower is None:
            f_lower = max(self.metadata.f_lower, waveform.metadata.f_lower)

        flen = 1 << (max(len(wf1), len(wf2)) - 1).bit_length()
        delta_f = 1.0 / (flen * delta_t)

        psd = aLIGOZeroDetHighPower(flen, delta_f, f_lower)
        wf1.resize(flen)
        wf2.resize(flen)

        return cbcmatch(wf1, wf2, psd=psd, low_frequency_cutoff=f_lower)[0]

    def freq(self, mode: tuple[int, int] = (2, 2)) -> np.ndarray:

        delta_t = self.metadata.delta_t

        wf_strain = 0
        for mode in self.metadata.modes:
            single_mode = self[mode]
            single_minus_mode = self[mode[0], -mode[1]]

            wf_strain += single_mode * spherical_harmonics(
                mode[0], mode[1], self.metadata.inclination, self.metadata.coa_phase
            ) + single_minus_mode * spherical_harmonics(
                mode[0], -mode[1], self.metadata.inclination, self.metadata.coa_phase
            )

        wf = ts.TimeSeries(wf_strain.real, delta_t=delta_t)

        wfreq = wf.to_frequencyseries()

        return wfreq

    def add_eccentricity(self, func, modes=[(2, 2)], **kwargs):
        strain = []
        for mode in modes:
            time, phase, amplitude = func(self, mode, **kwargs)
            strain.append(amplitude * np.exp(1j * phase))
        return BaseWaveform(np.vstack(strain), time, 0)
