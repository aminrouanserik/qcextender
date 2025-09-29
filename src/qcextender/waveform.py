import numpy as np
from typing import Iterable, Sequence, Self
from pycbc.waveform import get_td_waveform
from pycbc.psd import aLIGOZeroDetHighPower
from pycbc.types import timeseries as ts
from pycbc.filter.matchedfilter import match as cbcmatch
from qcextender.metadata import Metadata
from qcextender.basewaveform import BaseWaveform
from qcextender.utils import lal_waves


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
    ) -> Self:  # Dimensioned
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
                time, strain = lal_waves(
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

        # single_mode_strain = []
        # for mode in modes:
        #     local_kwargs = kwargs.copy()
        #     local_kwargs["mode_array"] = mode
        #     hp, hc = get_td_waveform(**local_kwargs)
        #     single_mode_strain.append((hp - 1j * hc))

        #     # single_mode_strain.append(
        #     #     waveform_modes.get_td_waveform_modes(**local_kwargs)
        #     # )

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
            float: The real match of the two waveforms.
        """

        wf1 = ts.TimeSeries(self[2, 2].real, delta_t=self.metadata.delta_t)
        wf2 = ts.TimeSeries(waveform[2, 2].real, delta_t=waveform.metadata.delta_t)

        if f_lower is None:
            f_lower = max(self.metadata.f_lower, waveform.metadata.f_lower)

        flen = 1 << (max(len(wf1), len(wf2)) - 1).bit_length()
        delta_f = 1.0 / (flen * wf1.delta_t)

        psd = aLIGOZeroDetHighPower(flen, delta_f, f_lower)
        wf1.resize(flen)
        wf2.resize(flen)

        return cbcmatch(wf1, wf2, psd=psd, low_frequency_cutoff=f_lower)[0]
