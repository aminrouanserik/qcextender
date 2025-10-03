from math import factorial, comb
import numpy as np
from numbers import Number


def safe_comb(n: int, k: int) -> int:
    return comb(n, k) if 0 <= k <= n else 0


def spherical_harmonics(l: int, m: int, iota: Number, phi: Number, s: int = -2):
    terms = [
        (-1) ** r
        * safe_comb(l - s, r)
        * safe_comb(l + s, r + s - m)
        * (1 / np.tan(iota / 2)) ** (2 * r + s - m)
        for r in range(0, l - s + 1)
    ]

    prefactor = (
        (-1) ** (l + m - s)
        * np.sqrt(
            (factorial(l + m) * factorial(l - m) * (2 * l + 1))
            / (4 * np.pi * factorial(l + s) * factorial(l - s))
        )
        * (np.sin(iota / 2)) ** (2 * l)
        * np.exp(1j * m * phi)
    )

    return prefactor * np.sum(terms)
