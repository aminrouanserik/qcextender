"""Module for loading, storing, and converting dimensionless gravitational waveforms.

This module provides the `DimensionlessWaveform` class, which wraps waveforms
from numerical relativity simulations (e.g., the SXS catalog). It stores
dimensionless time-domain data with associated metadata and allows conversion
to dimensional waveforms in SI units.

Conversions are handled using utility functions from `qcextender.units`, and
metadata is managed using the `Metadata` class. Dimensionless waveforms can be
converted to `Waveform` objects for direct use in physical analyses and model
comparisons.

Example:
    >>> from qcextender.dimensionlesswaveform import DimensionlessWaveform
    >>> dwf = DimensionlessWaveform.from_sim("SXS:BBH:1155")
    >>> wf = dwf.to_Waveform(f_lower=20, total_mass=60, distance=400)
"""

import sxs
import numpy as np
from typing import Self
from numbers import Number
from qcextender.metadata import Metadata
from qcextender.waveform import Waveform
from qcextender.basewaveform import BaseWaveform
from qcextender.units import tM_to_tSI, mM_to_mSI
from qcextender.functions import frequency_window, amp, phase


class DimensionlessWaveform(BaseWaveform):
    """Represents a dimensionless gravitational waveform from an NR simulation.

    The waveform contains multiple spherical harmonic modes and associated simulation metadata. It can be converted into a dimensional
    `Waveform` object using the `to_Waveform` method.

    Attributes:
        strain (np.ndarray): Stacked complex strain data for each mode.
        time (np.ndarray): Time array corresponding to the waveform.
        metadata (Metadata): Object storing waveform parameters and provenance.
    """

    def __init__(
        self, strain: np.ndarray, time: np.ndarray, metadata: Metadata
    ) -> None:
        """Initialize a dimensionless waveform.

        Args:
            strain (np.ndarray): Stacked array of complex strain modes.
            time (np.ndarray): Time array (dimensionless units), same length as strain arrays.
            metadata (Metadata): Simulation metadata describing the waveform.
        """
        super().__init__(strain, time, metadata)

    @classmethod
    def from_sim(cls, sim_id: str, modes: list[tuple[int, int]] = [(2, 2)]) -> Self:
        """Load a dimensionless waveform from an SXS simulation.

        Args:
            sim_id (str): Simulation identifier in the SXS catalog (e.g., "SXS:BBH:1155").
            modes (list[tuple[int, int]]): List of (l, m) modes to include. Defaults to [(2, 2)].

        Raises:
            ValueError: If a requested mode is not available in the simulation.

        Returns:
            DimensionlessWaveform: Waveform object containing the requested modes.
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
        """Convert the dimensionless waveform to a dimensional `Waveform`.

        The conversion scales time and strain using physical parameters such as total mass and distance. The resulting waveform is cropped to the longest
        continuous segment where the orbital frequency exceeds `f_lower`.

        Args:
            f_lower (Number): Lower frequency bound of the signal [Hz].
            total_mass (Number): Total mass of the binary [solar masses].
            distance (Number): Luminosity distance to the source [megaparsecs].
            inclination (Number, optional): Inclination angle of the system [radians]. Defaults to 0.
            coa_phase (Number, optional): Coalescence phase [radians]. Defaults to 0.

        Raises:
            ValueError: If no part of the waveform remains above `f_lower`.

        Returns:
            Waveform: Dimensional waveform object with physical units.
        """
        time = tM_to_tSI(self.time, total_mass)
        metadata = self.metadata.copy()
        newmetadata = metadata.to_dimensional(
            f_lower, total_mass, distance, inclination, coa_phase
        )

        single_mode_strains = []
        for mode in self.metadata.modes:
            singlemode = mM_to_mSI(self[mode], total_mass, distance)

            strain, time = frequency_window(singlemode, time, f_lower)

            phases, amps = phase(strain), amp(strain)
            phases -= phases[np.argmax(amps)]
            single_mode_strains.append(amps * np.exp(1j * phases))

        strain = np.vstack(single_mode_strains)
        return Waveform(strain, time, newmetadata)
