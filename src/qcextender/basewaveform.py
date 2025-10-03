import numpy as np
import spheroidal
from qcextender.metadata import Metadata
from numbers import Number
from dataclasses import fields


class BaseWaveform:
    """Base for all Waveform objects, contains all necessary attributes and methods."""

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
            index = modes.index((mode[0], abs(mode[1])))
        except ValueError:
            raise ValueError(f"Mode {mode} not found in this waveform.")

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
        return spheroidal.sphericalY(-2, l, m)(iota, phi)

    @staticmethod
    def _kwargs_to_metadata(
        kwargs: dict[str, type],
    ) -> "Metadata":
        """Converts SXS metadata or user supplied kwargs into a uniform Metadata object. Will be split for SXS and PyCBC later.

        Args:
            kwargs (dict[str, type]): Keyword arguments used to generate waveform or SXS simulation metadata.

        Returns:
            Metadata: Unified object encoding all important metadata.
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

    def abs(self, mode: tuple[int, int] = [2, 2]) -> np.ndarray:
        """Returns the amplitude of the mode strain.

        Args:
            mode (tuple[int, int], optional): Mode of which the amplitude is requested. Defaults to the dominant mode [2, 2].

        Returns:
            np.ndarray: The amplitude of the mode.
        """
        return np.abs(self[mode])

    def omega(self, mode: tuple[int, int] = [2, 2]) -> np.ndarray:
        """Returns the omega for a single mode.

        Args:
            mode (tuple[int, int], optional): Mode of which the omega is requested. Defaults to the dominant mode [2, 2].

        Returns:
            np.ndarray: The omega of the mode.
        """
        return np.gradient(-np.unwrap(np.angle(self[mode])), self.time)
