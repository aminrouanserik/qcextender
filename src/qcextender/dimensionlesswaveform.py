import sxs
import lal
import numpy as np
from typing import Iterable, Sequence, Self
from numbers import Number
from qcextender.metadata import Metadata
from qcextender.waveform import Waveform
from qcextender.basewaveform import BaseWaveform
from scipy.interpolate import make_interp_spline


class DimensionlessWaveform(BaseWaveform):
    """Dimensionless waveform class which can be converted to a standard Waveform class, stores multi-modal (time domain) waveforms."""

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

    @staticmethod
    def _strain_m(strain: np.ndarray, mass: Number, distance: Number) -> np.ndarray:
        """Converts geometric strain into strain in meters.

        Args:
            strain (np.ndarray): N-dimensional array containing strain in geometric units.
            mass (Number): Total mass in solar masses
            distance (Number): Distance in Mpc.

        Returns:
            np.ndarray: Wave strain in meters.
        """
        distance *= 1e6 * lal.PC_SI
        correction = mass * lal.MTSUN_SI * lal.C_SI / distance
        return strain * correction

    @staticmethod
    def _time_sec(time: np.ndarray, mass: Number) -> np.ndarray:
        """Converts geometric time into time in seconds.

        Args:
            time (np.ndarray): Time in geometric units.
            mass (Number): Mass in solar masses.

        Returns:
            np.ndarray: Time in seconds.
        """
        return time * (lal.MTSUN_SI * mass)

    @staticmethod
    def _freq_Hz(frequency, mass):
        return frequency / (lal.MTSUN_SI * mass)

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
        time = self._time_sec(self.time, total_mass)

        single_mode_strains = []
        for strain in self.strain:
            singlemode = self._strain_m(strain, total_mass, distance)

            # Scales decently now, probably not the most reliable and still see some issues
            freq = self._freq_Hz(
                np.gradient(-np.unwrap(np.angle(self.strain[0])), self.time), total_mass
            )
            arg = np.abs(singlemode)
            phase = np.unwrap(np.angle(singlemode))

            try:
                # Need a better way of doing this, too dependent on the one example I saw.
                # Why 4 * f_lower?
                cutoff = int(np.argwhere(np.isclose(freq, 4 * f_lower, atol=1))[-2])
            except:
                cutoff = 0

            interpolatedarg = make_interp_spline(time[cutoff:], arg[cutoff:])(
                time[cutoff:]
            )
            # Probably need something smarter than this, start is too sudden instead of gradual like in models
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
