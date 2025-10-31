"""Base waveform module for the qcextender package.

Defines the `BaseWaveform` class, which implements common data structures and
operations for all waveform objects. It provides utilities for:
- Mode access and recombination.
- Time alignment.
- Metadata normalization.
- Amplitude, phase, and frequency (omega) computation.

This class is subclassed by `Waveform` and `DimensionlessWaveform` to provide
specific domain behavior (e.g., dimensional vs. dimensionless representations).
"""

import numpy as np
from scipy.interpolate import make_interp_spline
from dataclasses import fields
from qcextender.metadata import Metadata
from qcextender.functions import spherical_harmonics, amp, phase, omega


class BaseWaveform:
    """Base for all Waveform objects, containing all core attributes and methods.

    This class handles generic waveform functionality such as time alignment, strain recombination, and metadata normalization. It is not meant to be
    instantiated directly but serves as the foundation for higher-level classes like `Waveform` and `DimensionlessWaveform`.

    Attributes:
        strain (np.ndarray): Stacked complex strain data for each mode.
        time (np.ndarray): Time array corresponding to the waveform.
        metadata (Metadata): Object storing waveform parameters and provenance.
    """

    def __init__(
        self, strain: np.ndarray, time: np.ndarray, metadata: Metadata
    ) -> None:
        """Initializes class containing waveform data.

        Args:
            strain (np.ndarray): Stacked array of multi-modal wave strains.
            time (np.ndarray): Time array, should be the same length as component strain arrays.
            metadata (Metadata): Metadata belonging to the generated or requested waveform.
        """
        self.strain = strain
        self.time = time
        self.metadata = metadata

    def __getitem__(self, mode: tuple[int, int]) -> np.ndarray:
        """Returns the single-mode wave strain corresponding to (l, m).

        Args:
            mode (tuple[int, int]): Spherical harmonics decomposed strain mode.

        Raises:
            ValueError: If this waveform object does not contain the requested mode.

        Returns:
            np.ndarray: The complex strain for the requested mode.
        """
        modes = self.metadata["modes"]
        try:
            index = modes.index((mode[0], abs(mode[1])))
        except ValueError:
            raise ValueError(f"Mode {mode} not found in this waveform.")

        # For negative m, return the conjugate mode with parity correction
        if mode[1] < 0:
            return (-1) ** mode[0] * np.conj(self.strain[index])
        return self.strain[index]

    @staticmethod
    def _align(strain: np.ndarray, time: np.ndarray) -> np.ndarray:
        """Aligns waveform such that the dominant order peak is at t=0.

        Args:
            strain (np.ndarray): Waveform strain data.
            time (np.ndarray): Time array prior to realigning.

        Returns:
            np.ndarray: New time array aligned to the maximum amplitude.
        """
        time -= time[np.argmax(np.abs(strain))]
        return time

    @staticmethod
    def _kwargs_to_metadata(kwargs: dict[str, type]) -> Metadata:
        """Converts simulation metadata or kwargs into a uniform Metadata object.

        Args:
            kwargs (dict[str, type]): Keyword arguments used to generate waveform or SXS metadata.

        Returns:
            Metadata: Unified object encoding all important metadata fields.
        """
        meta_fields = {f.name for f in fields(Metadata)}

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

    def recombine_strain(self, time: np.ndarray | None = None) -> np.ndarray:
        """Recombines individual (l, m) strain modes into a total observed strain.

        Args:
            time (np.ndarray | None, optional): Optional new time array for resampling via spline interpolation.
                If None, uses the waveform's native time array.

        Returns:
            np.ndarray: Complex strain representing the full waveform at the given inclination and phase.
        """
        strain = 0
        for mode in self.metadata.modes:
            if time is not None:
                single_mode = make_interp_spline(self.time, self[mode])(time)
                single_minus_mode = make_interp_spline(
                    self.time, self[mode[0], -mode[1]]
                )(time)
            else:
                single_mode = self[mode]
                single_minus_mode = self[mode[0], -mode[1]]

            strain += single_mode * spherical_harmonics(
                mode[0], mode[1], self.metadata.inclination, self.metadata.coa_phase
            ) + single_minus_mode * spherical_harmonics(
                mode[0],
                -mode[1],
                self.metadata.inclination,
                self.metadata.coa_phase,
            )
        return strain

    def amp(self, mode: tuple[int, int] = (2, 2)) -> np.ndarray:
        """Returns the amplitude of the specified mode.

        Args:
            mode (tuple[int, int], optional): Mode for which the amplitude is requested. Defaults to (2, 2).

        Returns:
            np.ndarray: The amplitude of the mode.
        """
        return amp(self[mode])

    def phase(self, mode: tuple[int, int] = (2, 2)) -> np.ndarray:
        """Returns the phase of the specified mode.

        Args:
            mode (tuple[int, int], optional): Mode for which the phase is requested. Defaults to (2, 2).

        Returns:
            np.ndarray: The phase of the mode.
        """
        return phase(self[mode])

    def omega(self, mode: tuple[int, int] = (2, 2)) -> np.ndarray:
        """Returns the angular frequency (omega) of the specified mode.

        Args:
            mode (tuple[int, int], optional): Mode for which the omega is requested. Defaults to (2, 2).

        Returns:
            np.ndarray: The instantaneous angular frequency of the mode.
        """
        return omega(self[mode], self.time)
