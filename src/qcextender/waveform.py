from typing import Iterable, Sequence, Self, Optional
from numbers import Number
import numpy as np
from pycbc.waveform import get_td_waveform, waveform_modes
from pycbc.types import timeseries as ts
from pycbc.filter.matchedfilter import match as cbcmatch
from pycbc.psd import aLIGOZeroDetHighPower
import sxs
import lal
from dataclasses import dataclass, field, fields, asdict


@dataclass
class Metadata:

    library: str
    q: float

    approximant: Optional[str] = None
    simulation_id: Optional[str] = None

    total_mass: Optional[float] = None
    spin1: tuple[float, float, float] = (0, 0, 0)
    spin2: tuple[float, float, float] = (0, 0, 0)
    eccentricity: float = 0
    distance: Optional[float] = None
    inclination: float = 0
    coa_phase: float = 0

    delta_t: Optional[float] = 1.0 / 4096
    f_lower: float = 20
    # f_final: Optional[float] = None

    modes: Iterable[tuple[int, int]] = field(default_factory=list)
    domain: str = "time"
    dimensionless: bool = True
    aligned_to_peak: bool = True

    def __getitem__(self, key):
        return getattr(self, key)

    def to_dict(self):
        return asdict(self)


class Waveform:
    """Base waveform class with which all calculations can be done, stores multi-modal (time domain) waveforms."""

    def __init__(
        self, strain: np.ndarray, time: np.ndarray, metadata: Metadata
    ) -> None:
        """Initializes class containing waveform data.

        Args:
            strain (np.ndarray): Stacked array of multi-modal wave strains.
            time (np.ndarray): Time array, should be the same length as component strain arrays.
            metadata (dict): Metadata belonging to the generated or requested waveform.
        """
        self.strain = strain
        self.time = time
        self.metadata = metadata

    @classmethod
    def from_model(
        cls, approximant: str, modes: Iterable[Sequence[int]], **kwargs
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
            library="PyCBC",
            q=q,
            approximant=approximant,
            modes=list(modes),
        )

        single_mode_strain = []
        for mode in modes:
            local_kwargs = kwargs.copy()
            local_kwargs["mode_array"] = mode
            hp, hc = get_td_waveform(**local_kwargs)
            single_mode_strain.append(np.asarray(hp - 1j * hc))

        multi_mode_strain = cls._dimensionless_strain(
            np.vstack(single_mode_strain),
            total_mass,
            kwargs["distance"],
        )
        kwargs["distance"] = None
        time = cls._dimensionless_time(hp.sample_times, total_mass)
        metadata = cls._kwargs_to_metadata(kwargs)

        return cls(
            multi_mode_strain,
            time,
            metadata,
        )

    @classmethod
    def from_sim(cls, sim_id: str, modes: Iterable[Sequence[int]]) -> Self:
        """Loads multi-modal data from a specified SXS simulation.

        Args:
            name (str): Name of the simulation as in the SXS catalogue.
            modes (Iterable[Sequence[int]]): Modes to be included.

        Raises:
            ValueError: A specified mode in modes is not available.

        Returns:
            Self: Object with stacked modes as complex array.
        """
        sim = sxs.load(sim_id, extrapolation="Outer")
        meta = sim.metadata
        meta["modes"] = modes

        q = meta["initial_mass1"] / meta["initial_mass2"]
        if q < 1:
            q = 1 / q

        meta.update(
            library="SXS",
            simulation_id=sim_id,
            q=q,
            modes=list(modes),
        )

        single_mode_strain = []
        for l, m in modes:
            try:
                single_mode_strain.append(
                    np.array(sim.h[:, sim.h.index(l, m)])
                    * cls._spherical_harmonic(
                        l, m, 0, np.pi / 2
                    )  # Shifted to make visual comparisons easier
                )
            except:
                raise ValueError(f"Mode (l={l}, m={m}) not found in this simulation.")

        time = cls._align(np.array(sim.h[:, sim.h.index(l, m)]), sim.h.t)
        multi_mode_strain = np.vstack(single_mode_strain)
        metadata = cls._kwargs_to_metadata(meta)
        return cls(multi_mode_strain, time, metadata)

    @staticmethod
    def _dimensionless_strain(
        strain: np.ndarray, mass: Number, distance: Number
    ) -> np.ndarray:
        """Converts the strain into a dimensionless strain.

        Args:
            strain (np.ndarray): N-dimensional array containing dimensioned strains.
            mass (Number): Total mass in solar masses
            distance (Number): Distance in Mpc.

        Returns:
            np.ndarray: Dimensionless wave strain.
        """
        distance *= 1e6 * lal.PC_SI
        correction = mass * lal.MTSUN_SI * lal.C_SI / distance
        return strain / correction

    @staticmethod
    def _dimensionless_time(time: np.ndarray, mass: Number) -> np.ndarray:
        """Converts time from units of seconds to geometric units.

        Args:
            time (np.ndarray): Time in seconds.
            mass (Number): Mass in solar masses.

        Returns:
            np.ndarray: New time in geometric units.
        """
        return time / (lal.MTSUN_SI * mass)

    @staticmethod
    def _align(strain: np.ndarray, time: np.ndarray) -> np.ndarray:
        """Aligns waveform such that the dominant order peak is at t=0.

        Args:
            strain (np.ndarray): Waveform strain data.
            time (np.ndarray): Time array prior to realigning.

        Returns:
            np.ndarray: New time array.
        """
        time -= time[np.argmax(np.abs(strain))]
        return time

    @staticmethod
    def _spherical_harmonic(l: Number, m: Number, iota: Number, phi: Number) -> Number:
        """Calculates the spin-weighted spherical harmonics for any mode.

        Args:
            l (Number): Spherical harmonic degree.
            m (Number): Spherical harmonic order.
            iota (Number): The inclination in radians.
            phi (Number): The coalescence phase.

        Returns:
            Number: The spin-weighted spherical harmonics at specified order.
        """
        return waveform_modes.get_glm(l, m, iota) * np.exp(1j * m * phi)

    @staticmethod
    def _kwargs_to_metadata(kwargs: dict[str, type]) -> "Metadata":
        """Converts SXS metadata or user supplied kwargs into a uniform Metadata object. Will be split for SXS and PyCBC later.

        Args:
            kwargs (dict[str, type]): Keyword arguments used to generate waveform or SXS simulation metadata.

        Returns:
            Metadata: Unified object encoding all important metadata.
        """
        meta_fields = {f.name for f in fields(Metadata)}

        kwargs["total_mass"] = None
        kwargs["distance"] = None

        aliases = {
            "reference_dimensionless_spin1": "spin1",
            "reference_dimensionless_spin2": "spin2",
            "reference_eccentricity": "eccentricity",
        }

        fixed_kwargs = {}
        for k, v in kwargs.items():
            k = aliases.get(k, k)
            if k in meta_fields:
                if k in ("spin1", "spin2") and v is not None:
                    fixed_kwargs[k] = tuple(v)  # normalize to tuple
                else:
                    fixed_kwargs[k] = v

        return Metadata(**fixed_kwargs)

    def __getitem__(self, mode: tuple[int, int]) -> np.ndarray:
        """Returns the single mode wave strain.

        Args:
            mode (tuple[int, int]): Spherical harmonics decomposed strain mode.

        Raises:
            ValueError: This Waveform object does not contain the requested mode.

        Returns:
            np.ndarray: Single mode wave strain.
        """

        modes = self.metadata["modes"]
        try:
            index = modes.index((mode[0], mode[1]))
        except ValueError:
            raise ValueError(f"Mode {mode} not found in this waveform.")
        return self.strain[index]

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
