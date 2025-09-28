from dataclasses import dataclass, field, asdict, fields
from typing import Iterable, Optional


@dataclass
class Metadata:
    """Class that unifies metadata or kwargs of waveforms from different origins.

    Raises:
        ValueError: If initialized incorrectly because of an unexpected combination of values.

    Returns:
        Metadata: Metadata of the waveform.
    """

    library: str
    q: float
    delta_t: float

    approximant: Optional[str] = None
    simulation_id: Optional[str] = None

    total_mass: Optional[float] = None
    distance: Optional[float] = None

    spin1: tuple[float, float, float] = (0, 0, 0)
    spin2: tuple[float, float, float] = (0, 0, 0)
    eccentricity: float = 0
    inclination: float = 0
    coa_phase: float = 0

    f_lower: Optional[float] = None
    # f_final: Optional[float] = None

    modes: Iterable[tuple[int, int]] = field(default_factory=list)
    # domain: str = "time"
    dimensionless: bool = True
    aligned_to_peak: bool = True

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
        "delta_t",
        "modes",
        "dimensionless",
        "aligned_to_peak",
    ]

    def __getitem__(self, key: str) -> type:
        """Get an item from the metadata.

        Args:
            key (str): Name of the requested item.

        Returns:
            type: Current value corresponding to the key.
        """
        return getattr(self, key)

    def __repr__(self) -> str:
        """Representation of the metadata, only prints assigned values.

        Returns:
            str: Representation of the Metadata object as a string.
        """
        parts = []
        for name in self._REPR_ORDER:
            val = getattr(self, name, None)
            if val is not None:
                parts.append(f"{name}={val}")
        return f"Metadata({', '.join(parts)})"

    def __post_init__(self):
        """Runs after initialization enforcing proper (non) dimensionality.

        Raises:
            ValueError: If missing keys that are crucial to dimensionality.
        """
        if self.dimensionless:
            self.total_mass = None
            self.distance = None
        else:
            if self.total_mass == None or self.distance == None:
                raise ValueError(
                    "Both total_mass and distance have to be defined in the dimensioned case."
                )

    def to_dict(self) -> dict:
        """Returns Metadata as a dictionary.

        Returns:
            dict: Metadata as a dictionary.
        """
        return asdict(self)
