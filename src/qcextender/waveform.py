from typing import Iterable, Sequence, Self
import numpy as np
from pycbc.waveform import get_td_waveform


class Waveform:
    """Base waveform class with which all calculations can be done, stores multi-modal (time domain) waves."""

    def __init__(self, strain: np.ndarray, time: np.ndarray, metadata: dict) -> None:
        """Agnostic class initializer, can supplement with any additional model."""
        self.strain = strain
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
            single_mode_strain.append(np.asarray(hp + 1j * hc))

        multi_mode_strain = np.vstack(single_mode_strain)
        return cls(multi_mode_strain, hp.sample_times, kwargs)

    def singlemode(self, l: int = 2, m: int = 2) -> np.ndarray:
        """Returns the single mode wave strain. Defaults to the dominant mode.

        Args:
            l (int, optional): Spherical harmonic degree. Defaults to 2.
            m (int, optional): Spherical harmonic order. Defaults to 2.

        Raises:
            ValueError: The requested mode could not be found.

        Returns:
            np.ndarray: Single mode wave strain.
        """
        modes = self.metadata.get("modes", [])
        try:
            index = modes.index((l, m))
        except ValueError:
            raise ValueError(f"Mode (l={l}, m={m}) not found in this waveform.")
        return self.strain[index]
