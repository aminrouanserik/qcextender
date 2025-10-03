"""Utility module for converting between geometric units and SI units.

All functions assume geometric units with G = c = 1 unless otherwise noted.
"""

import numpy as np
from numbers import Number

MTSUN_SI = 4.925490947641267e-06
PC_SI = 3.085677581491367e16
C_SI = 299792458.0


def tM_to_tSI(time: np.ndarray | Number, total_mass: Number) -> np.ndarray:
    """Convert geometric time (units of M) into SI seconds.

    Args:
        time (array-like or float): Time in geometric units (M).
        total_mass (float): Total mass of the system [solar masses].

    Returns:
        np.ndarray: Time in seconds.
    """
    return time * (MTSUN_SI * total_mass)


def tSI_to_tM(time: np.ndarray | Number, total_mass: Number) -> np.ndarray:
    """Convert SI seconds into geometric time (units of M).

    Args:
        time (array-like or float): Time in seconds.
        total_mass (float): Total mass of the system [solar masses].

    Returns:
        np.ndarray: Time in geometric units (M).
    """
    return time / (MTSUN_SI * total_mass)


def fM_to_fSI(frequency: np.ndarray | Number, total_mass: Number) -> np.ndarray:
    """Convert geometric frequency (units of 1/M) into SI hertz.

    Args:
        frequency (array-like or float): Frequency in geometric units (1/M).
        total_mass (float): Total mass of the system [solar masses].

    Returns:
        np.ndarray: Frequency in hertz.
    """
    return frequency / (MTSUN_SI * total_mass)


def mM_to_mSI(
    strain: np.ndarray | Number, total_mass: Number, distance: Number
) -> np.ndarray:
    """Convert geometric strain amplitude into SI strain at a given distance.

    Args:
        strain (array-like or float): Dimensionless strain in geometric units.
        total_mass (float): Total mass of the system [solar masses].
        distance (float): Luminosity distance to source [megaparsecs].

    Returns:
        np.ndarray: Strain in SI units (meters).
    """
    return strain * total_mass * MTSUN_SI * C_SI / (distance * 1e6 * PC_SI)
