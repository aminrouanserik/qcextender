from typing import Iterable, Sequence, Self
from numbers import Number
import numpy as np
from pycbc.waveform import get_td_waveform, waveform_modes
import sxs
import lal


class Waveform:
    """Base waveform class with which all calculations can be done, stores multi-modal (time domain) waveforms."""

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

    def __getitem__(self, mode: tuple[int, int]) -> np.ndarray:
        """Returns the single mode wave strain.

        Args:
            mode (tuple[int, int]): Spherical harmonics decomposed strain mode.

        Raises:
            ValueError: This Waveform object does not contain the requested mode.

        Returns:
            np.ndarray: Single mode wave strain.
        """

        modes = self.metadata.get("modes", [])
        try:
            index = modes.index((mode[0], mode[1]))
        except ValueError:
            raise ValueError(f"Mode {mode} not found in this waveform.")
        return self.strain[index]
