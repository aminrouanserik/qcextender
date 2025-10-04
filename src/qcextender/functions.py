"""Core mathematical functions for gravitational waveform analysis.

Provides basic utilities to extract amplitude, phase, and instantaneous
frequency from complex gravitational-wave strains. These functions are
lightweight wrappers around NumPy operations and are intended for use
throughout the `qcextender` package.

Example:
    >>> from qcextender import functions as fn
    >>> strain = np.exp(1j * np.linspace(0, 4*np.pi, 1000))
    >>> fn.amp(strain)[0], fn.phase(strain)[-1]
    (1.0, 12.566370614359172)
"""

from math import factorial
import numpy as np
from scipy.special import binom
from numbers import Number


def amp(strain: np.ndarray) -> np.ndarray:
    """Compute the amplitude of a complex strain signal.

    Args:
        strain (np.ndarray): Complex gravitational-wave strain.

    Returns:
        np.ndarray: Amplitude |h(t)|.
    """
    return np.abs(strain)


def phase(strain: np.ndarray) -> np.ndarray:
    """Compute the unwrapped phase of a complex strain signal.

    Args:
        strain (np.ndarray): Complex gravitational-wave strain.

    Returns:
        np.ndarray: Unwrapped phase φ(t) in radians.
    """
    return np.unwrap(np.angle(strain))


def omega(strain: np.ndarray, time: np.ndarray) -> np.ndarray:
    """Compute the instantaneous angular frequency of a strain signal.

    Calculated as the time derivative of the unwrapped phase.

    Args:
        strain (np.ndarray): Complex gravitational-wave strain.
        time (np.ndarray): Time array corresponding to the strain samples [s].

    Returns:
        np.ndarray: Instantaneous angular frequency ω(t) [rad/s].
    """
    return np.gradient(-np.unwrap(np.angle(strain)), time)


def spherical_harmonics(
    l: int, m: int, iota: Number, phi: Number, s: Number = -2
) -> complex:
    """Calculates the spin-weighted spherical harmonics. Code adapted from the Spheroidal package.

    Args:
        l (int): Degree of the spherical harmonics.
        m (int): Integer order of the spherical harmonics.
        iota (Number): The inclination angle in radians.
        phi (Number): The phase of coalescence in radians.
        s (Number, optional): Int or half-integer spin-weight. Defaults to -2.

    Returns:
        complex: The spin-weighted spherical harmonics of mode (l, m)
    """
    prefactor = (-1.0) ** (l + m - s + 0j)
    prefactor *= np.sqrt(
        factorial(l + m)
        * factorial(l - m)
        * (2 * l + 1)
        / (4 * np.pi * factorial(l + s) * factorial(l - s))
    )

    alternating_sum = 0
    for r in range(int(max(m - s, 0)), int(min(l - s, l + m) + 1)):
        alternating_sum += (
            (-1) ** r
            * binom(l - s, r)
            * binom(l + s, r + s - m)
            * np.sin(iota / 2) ** (2 * l - 2 * r - s + m)
            * np.cos(iota / 2) ** (2 * r + s - m)
        )

    return prefactor * np.exp(1j * m * phi) * alternating_sum
