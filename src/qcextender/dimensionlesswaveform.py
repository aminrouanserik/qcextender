import sxs
import numpy as np
from typing import Iterable, Sequence, Self
from numbers import Number
from qcextender.metadata import Metadata
from qcextender.waveform import Waveform
from qcextender.basewaveform import BaseWaveform
from qcextender.units import tM_to_tSI, mM_to_mSI
from scipy.interpolate import make_interp_spline


class DimensionlessWaveform(BaseWaveform):
    """Dimensionless waveform class which can be converted to a standard Waveform class, stores multiple modes, time, and simulation metadata."""

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
    def from_sim(cls, sim_id: str, modes: Iterable[Sequence[int]] = [(2, 2)]) -> Self:
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

        dt = np.min(np.diff(sim.h.t))
        sim = sim.h.interpolate(np.arange(sim.h.t[0], sim.h.t[-1], dt))

        meta.update(
            library="SXS",
            simulation_id=sim_id,
            q=q,
            modes=list(modes),
            delta_t=dt,
            dimensionless=True,
        )

        single_mode_strain = []
        for l, m in modes:
            try:
                single_mode_strain.append(np.array(sim[:, sim.index(l, m)]))
            except:
                raise ValueError(f"Mode (l={l}, m={m}) not found in this simulation.")

        time = cls._align(np.array(sim[:, sim.index(l, m)]), sim.t)
        multi_mode_strain = np.vstack(single_mode_strain)
        metadata = cls._kwargs_to_metadata(meta)
        return cls(multi_mode_strain, time, metadata)

    def to_Waveform(
        self,
        f_lower: Number,
        total_mass: Number,
        distance: Number,
        inclination: Number = 0,
        coa_phase: Number = 0,
    ) -> Waveform:
        """Creates a copy of any dimensionless waveform and casts to a regular waveform with the specified parameters.

        Args:
            f_lower (Number): Lower frequency bound of the signal.
            total_mass (Number): Total mass of the binary in solar mass.
            distance (Number): Distance to the merger in Mpc.
            inclination (Number, optional): Inclination angle of the inspiral. Defaults to 0.
            coa_phase (Number, optional): Coalescence phase. Defaults to 0.

        Returns:
            Waveform: Waveform object with the admitted properties.
        """
        time = tM_to_tSI(self.time, total_mass)

        single_mode_strains = []
        for mode in self.metadata.modes:
            singlemode = mM_to_mSI(self[mode], total_mass, distance)

            omega = np.gradient(-np.unwrap(np.angle(singlemode)), time)
            arg = np.abs(singlemode)
            phase = np.unwrap(np.angle(singlemode))

            # Index of the argwhere needs to be redone and tested, reasonable now
            try:
                fpoints = np.argwhere(np.isclose(omega, 2 * np.pi * f_lower, atol=0.1))
                cutoff = int(fpoints[len(fpoints) // 2])
            except:
                cutoff = 0

            interpolatedarg = make_interp_spline(time[cutoff:], arg[cutoff:])(
                time[cutoff:]
            )
            interpolatedphase = make_interp_spline(time[cutoff:], phase[cutoff:])(
                time[cutoff:]
            )
            single_mode_strains.append(interpolatedarg * np.exp(1j * interpolatedphase))

        strain = np.vstack(single_mode_strains)
        metadata = self.metadata.copy()
        newmetadata = metadata.to_dimensional(
            f_lower, total_mass, distance, inclination, coa_phase
        )

        return Waveform(strain, time[cutoff:], newmetadata)
