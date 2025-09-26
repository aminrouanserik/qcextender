from typing import Iterable, Sequence, Self
from numbers import Number
import numpy as np
from pycbc.waveform import get_td_waveform
import sxs
import lal


class Waveform:
    """Base waveform class with which all calculations can be done, stores multi-modal (time domain) waves."""

    def __init__(self, strain: np.ndarray, time: np.ndarray, metadata: dict) -> None:
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
        single_mode_strain = []
        kwargs["approximant"] = approximant
        kwargs["modes"] = modes

        for mode in modes:
            local_kwargs = kwargs.copy()
            local_kwargs["mode_array"] = mode
            hp, hc = get_td_waveform(**local_kwargs)
            single_mode_strain.append(np.asarray(hp - 1j * hc))

        multi_mode_strain = cls._dimensionless_strain(
            np.vstack(single_mode_strain),
            kwargs["mass1"] + kwargs["mass2"],
            kwargs["distance"],
        )
        time = cls._dimensionless_time(
            hp.sample_times, kwargs["mass1"] + kwargs["mass2"]
        )

        return cls(
            multi_mode_strain,
            time,
            kwargs,
        )

    @classmethod
    def from_sim(cls, name: str, modes: Iterable[Sequence[int]]) -> Self:
        """Loads multi-modal data from a specified SXS simulation.

        Args:
            name (str): Name of the simulation as in the SXS catalogue.
            modes (Iterable[Sequence[int]]): Modes to be included.

        Raises:
            ValueError: A specified mode in modes is not available.

        Returns:
            Self: Object with stacked modes as complex array.
        """
        sim = sxs.load(name, extrapolation="Outer")
        metadata = sim.metadata
        metadata["modes"] = modes

        single_mode_strain = []
        for l, m in modes:
            try:
                single_mode_strain.append(np.array(sim.h[:, sim.h.index(l, m)]))
            except:
                raise ValueError(f"Mode (l={l}, m={m}) not found in this simulation.")

        time = cls._align(np.array(sim.h[:, sim.h.index(l, m)]), sim.h.t)
        multi_mode_strain = np.vstack(single_mode_strain)

        return cls(multi_mode_strain, time, metadata)

    # Spherical harmonics should be moved to a function that takes l, m, inclination and coa_phase.
    @staticmethod
    def _dimensionless_strain(
        strain: np.ndarray, mass: Number, distance: Number
    ) -> np.ndarray:
        """Converts the strain into a dimensionless strain. Currently includes spherical harmonics.

        Args:
            strain (np.ndarray): N-dimensional array containing dimensioned strains.
            mass (Number): Total mass in solar masses
            distance (Number): Distance in Mpc.

        Returns:
            np.ndarray: Dimensionless wave strain.
        """
        distance *= 1e6 * lal.PC_SI
        y22 = np.sqrt(5.0 / (64 * np.pi)) * ((1 + np.cos(0)) ** 2) * np.exp(2 * 0 * 1j)
        correction = mass * lal.MTSUN_SI * lal.C_SI / distance * y22
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
        time -= time[np.argmax(strain)]
        return time

    def singlemode(self, l: int = 2, m: int = 2) -> np.ndarray:
        """Returns the single mode wave strain. Defaults to the dominant mode.

        Args:
            l (int, optional): Spherical harmonic degree. Defaults to 2.
            m (int, optional): Spherical harmonic order. Defaults to 2.

        Raises:
            ValueError: The requested mode could not be found.

        Returns:
            np.ndarray: Single mode wave strain.
        """
        modes = self.metadata.get("modes", [])
        try:
            index = modes.index((l, m))
        except ValueError:
            raise ValueError(f"Mode (l={l}, m={m}) not found in this waveform.")
        return self.strain[index]
