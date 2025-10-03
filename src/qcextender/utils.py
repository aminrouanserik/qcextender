from math import factorial
import numpy as np
from scipy.special import binom
from numbers import Number


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
        alternating_sum = alternating_sum + (-1) ** r * binom(l - s, r) * binom(
            l + s, r + s - m
        ) * np.sin(iota / 2) ** (2 * l - 2 * r - s + m) * np.cos(iota / 2) ** (
            2 * r + s - m
        )

    return prefactor * np.exp(1j * m * phi) * alternating_sum
