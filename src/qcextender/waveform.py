from typing import Iterable, Sequence, Self
import numpy as np
from pycbc.waveform import get_td_waveform


class Waveform:
    """Base waveform class with which all calculations can be done, stores multi-modal (time domain) waves."""

    def __init__(self, data: np.ndarray, time: np.ndarray, metadata: dict) -> None:
        """Agnostic class initializer, can supplement with any additional model."""
        self.data = data
        self.time = time
        self.metadata = metadata

    @classmethod
    def from_sim(cls, name: str, modes: Iterable[Sequence[int]]) -> None:
        """Loads multi-modal data from a simulation."""
        pass

    @classmethod
    def from_model(
        cls, approximant: str, modes: Iterable[Sequence[int]], **kwargs
    ) -> Self:
        """Generates a multi-modal time-domain waveform from a specified model using PyCBC.

        Args:
            approximant (str): PyCBC included waveform models.
            modes (Iterable[Sequence[int]]): Modes to be included.

        Returns:
            Self: _description_
        """
        kwargs["approximant"], kwargs["mode_array"] = approximant, modes
        hp, hc = get_td_waveform(**kwargs)
        return cls(np.asarray(hp + 1j * hc), hp.sample_times, kwargs)
