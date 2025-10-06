"""Metadata handling for waveform objects in the qcextender package.

This module defines the `Metadata` dataclass, a unified container for waveform
metadata originating from different sources (e.g., SXS simulations or lalsimulation models).

The Metadata class standardizes key waveform properties (mass ratio, spins,
frequency bounds, etc.) and provides methods for copying, converting between
dimensionless and dimensional representations, and dictionary serialization.
"""

from dataclasses import dataclass, field, asdict
from typing import Iterable, Optional, Self
from numbers import Number
from qcextender.units import tM_to_tSI


@dataclass
class Metadata:
    """Unified container for waveform metadata or generation parameters.

    This class ensures consistent handling of waveform metadata between
    simulation-based (dimensionless) and physical (dimensional) waveforms.

    Raises:
        ValueError: If initialized incorrectly due to missing required dimensional parameters.
    """

    # Core properties
    library: str
    q: float
    delta_t: float

    # Model/simulation identification
    approximant: Optional[str] = None
    simulation_id: Optional[str] = None

    # Physical parameters
    total_mass: Optional[float] = None
    distance: Optional[float] = None

    # Spin and orbital parameters
    spin1: tuple[float, float, float] = (0, 0, 0)
    spin2: tuple[float, float, float] = (0, 0, 0)
    eccentricity: float = 0
    inclination: float = 0
    coa_phase: float = 0

    # Frequency-related parameters
    f_lower: Optional[float] = None
    f_ref: Optional[float] = None

    # Mode and domain information
    modes: Iterable[tuple[int, int]] = field(default_factory=list)
    dimensionless: bool = False
    aligned_to_peak: bool = True

    # Preferred order for __repr__
    _REPR_ORDER = [
        "library",
        "approximant",
        "simulation_id",
        "q",
        "total_mass",
        "distance",
        "spin1",
        "spin2",
        "eccentricity",
        "inclination",
        "coa_phase",
        "f_lower",
        "f_ref",
        "delta_t",
        "modes",
        "dimensionless",
        "aligned_to_peak",
    ]

    def __getitem__(self, key: str) -> type:
        """Access a metadata attribute using dict-like indexing.

        Args:
            key (str): Name of the requested item.

        Returns:
            type: Value of the key.
        """
        return getattr(self, key)

    def __repr__(self) -> str:
        """String representation of the metadata.

        Only includes assigned values in a fixed order for readability.

        Returns:
            str: Formatted string representation of the Metadata object.
        """
        parts = []
        for name in self._REPR_ORDER:
            val = getattr(self, name, None)
            if val is not None:
                parts.append(f"{name}={val}")
        return f"Metadata({', '.join(parts)})"

    def __post_init__(self):
        """Post-initialization validation.

        Ensures correct consistency between dimensional and dimensionless metadata.

        Raises:
            ValueError: If required parameters are missing for dimensional metadata.
        """
        if self.dimensionless:
            self.total_mass = None
            self.distance = None
        else:
            if self.total_mass is None or self.distance is None:
                raise ValueError(
                    "Both total_mass and distance must be defined for dimensional waveforms."
                )

    def __copy__(self) -> Self:
        """Creates a shallow copy of this Metadata object.

        Returns:
            Self: Shallow copy of the Metadata instance.
        """
        cls = self.__class__
        copied = cls(
            **{
                f.name: getattr(self, f.name)
                for f in self.__dataclass_fields__.values()
            }
        )

        # Copy mutable attributes explicitly
        if isinstance(self.modes, list):
            copied.modes = self.modes.copy()

        return copied

    def copy(self) -> Self:
        """Alias for __copy__ to follow common Python convention.

        Returns:
            Self: Shallow copy of Metadata object.
        """
        return self.__copy__()

    def to_dict(self) -> dict:
        """Returns metadata as a dictionary.

        Returns:
            dict: Dictionary representation of this Metadata object.
        """
        return asdict(self)

    def to_dimensional(
        self,
        f_lower: Number,
        total_mass: Number,
        distance: Number,
        inclination: Number,
        coa_phase: Number,
    ) -> Self:
        """Converts this dimensionless Metadata instance into a dimensional one.

        This method updates frequency, mass, distance, and angle parameters, converts delta_t from M-units to SI seconds, and toggles the
        `dimensionless` flag to False. Copy metadata first to avoid overwriting.

        Args:
            f_lower (Number): Lower frequency bound of the waveform.
            total_mass (Number): Total mass of the system in solar masses.
            distance (Number): Distance to the source in Mpc.
            inclination (Number): Inclination angle of the system.
            coa_phase (Number): Coalescence phase.

        Returns:
            Self: The same Metadata instance, modified in-place and dimensionalized.
        """
        self.f_lower = f_lower
        self.total_mass = total_mass
        self.distance = distance
        self.inclination = inclination
        self.coa_phase = coa_phase
        self.delta_t = tM_to_tSI(self.delta_t, total_mass)
        self.dimensionless = False
        return self
