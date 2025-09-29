import sxs
import lal
import numpy as np
from typing import Iterable, Sequence, Self
from numbers import Number
from qcextender.metadata import Metadata
from qcextender.basewaveform import BaseWaveform


class DimensionlessWaveform(BaseWaveform):
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
        super().__init__(strain, time, metadata)

    @classmethod
    def from_sim(
        cls, sim_id: str, modes: Iterable[Sequence[int]] = [(2, 2)]
    ) -> Self:  # Dimensionless
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

        q = meta["initial_mass_ratio"]
        if q < 1:
            q = 1 / q

        meta.update(
            library="SXS",
            simulation_id=sim_id,
            q=q,
            modes=list(modes),
            delta_t=sim.h.t[1] - sim.h.t[0],
            dimensionless=True,
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
    ) -> np.ndarray:  # Conversion -> dimensionless, rewrite
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
    def _dimensionless_time(
        time: np.ndarray, mass: Number
    ) -> np.ndarray:  # Conversion -> dimensionless, rewrite
        """Converts time from units of seconds to geometric units.

        Args:
            time (np.ndarray): Time in seconds.
            mass (Number): Mass in solar masses.

        Returns:
            np.ndarray: New time in geometric units.
        """
        return time / (lal.MTSUN_SI * mass)
