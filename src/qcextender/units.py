"""Utility module for converting between geometric units and SI units.

All functions assume geometric units with G = c = 1 unless otherwise noted.
Conversions are provided for time, frequency, and strain. Each conversion
has a corresponding inverse, ensuring round-trip consistency.

Constants:
    MTSUN_SI (float): Solar mass in seconds.
    PC_SI (float): Parsec in meters.
    C_SI (float): Speed of light in m/s.

Example:
    >>> from qcextender import units
    >>> units.tM_to_tSI(100, 30)
    0.0147764728429238
"""

import numpy as np
from numbers import Number

MTSUN_SI: float = 4.925490947641267e-06
PC_SI: float = 3.085677581491367e16
C_SI: float = 299792458.0


def tM_to_tSI(time: np.ndarray | Number, total_mass: Number) -> np.ndarray:
    """Convert geometric time (units of M) into SI seconds.

    Inverse of `tSI_to_tM`.

    Args:
        time (array-like or float): Time in geometric units (M).
        total_mass (float): Total mass of the system [solar masses].

    Returns:
        np.ndarray: Time in seconds.
    """
    return time * (MTSUN_SI * total_mass)


def tSI_to_tM(time: np.ndarray | Number, total_mass: Number) -> np.ndarray:
    """Convert SI seconds into geometric time (units of M).

    Inverse of `tM_to_tSI`.

    Args:
        time (array-like or float): Time in seconds.
        total_mass (float): Total mass of the system [solar masses].

    Returns:
        np.ndarray: Time in geometric units (M).
    """
    return time / (MTSUN_SI * total_mass)


def fM_to_fSI(frequency: np.ndarray | Number, total_mass: Number) -> np.ndarray:
    """Convert geometric frequency (units of 1/M) into SI hertz.

    Inverse of `fSI_to_fM`.

    Args:
        frequency (array-like or float): Frequency in geometric units (1/M).
        total_mass (float): Total mass of the system [solar masses].

    Returns:
        np.ndarray: Frequency in hertz.
    """
    return frequency / (MTSUN_SI * total_mass)


def fSI_to_fM(frequency: np.ndarray | Number, total_mass: Number) -> np.ndarray:
    """Convert SI hertz into geometric frequency (units of 1/M).

    Inverse of `fM_to_fSI`.

    Args:
        frequency (array-like or float): Frequency in hertz.
        total_mass (float): Total mass of the system [solar masses].

    Returns:
        np.ndarray: Frequency in geometric units (1/M).
    """
    return frequency * (MTSUN_SI * total_mass)


def mM_to_mSI(
    strain: np.ndarray | Number, total_mass: Number, distance: Number
) -> np.ndarray:
    """Convert geometric strain amplitude into SI strain at a given distance.

    Inverse of `mSI_to_mM`.

    Args:
        strain (array-like or float): Dimensionless strain in geometric units.
        total_mass (float): Total mass of the system [solar masses].
        distance (float): Luminosity distance to source [megaparsecs].

    Returns:
        np.ndarray: Equivalent displacement amplitude in SI units (meters).
    """
    return strain * total_mass * MTSUN_SI * C_SI / (distance * 1e6 * PC_SI)


def mSI_to_mM(
    strain: np.ndarray | Number, total_mass: Number, distance: Number
) -> np.ndarray:
    """Convert SI strain amplitude into geometric strain at a given distance.

    Inverse of `mM_to_mSI`.

    Args:
        strain (array-like or float): Equivalent displacement amplitude in SI units (meters).
        total_mass (float): Total mass of the system [solar masses].
        distance (float): Luminosity distance to source [megaparsecs].

    Returns:
        np.ndarray: Dimensionless strain in geometric units.
    """
    return strain / (total_mass * MTSUN_SI * C_SI / (distance * 1e6 * PC_SI))
