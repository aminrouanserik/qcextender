"""Module for unit conversion functions as a utility."""

import numpy as np
from numbers import Number

MTSUN_SI = 4.925490947641267e-06
PC_SI = 3.085677581491367e16
C_SI = 299792458.0


def tM_to_tSI(time: np.ndarray, total_mass: Number) -> np.ndarray:
    """Converts geometric time into time in seconds.

    Args:
        time (np.ndarray): Geometric time.
        total_mass (Number): Total mass in solar masses.

    Returns:
        np.ndarray: Time in seconds.
    """
    return time * (MTSUN_SI * total_mass)


def fM_to_fSI(frequency: np.ndarray, total_mass: Number) -> np.ndarray:
    """Converts geometric frequency into frequency in hertz.

    Args:
        frequency (np.ndarray): Geometric frequency.
        total_mass (Number): Total mass in solar masses.

    Returns:
        np.ndarray: Frequency in hertz.
    """
    return frequency / (MTSUN_SI * total_mass)


def mM_to_mSI(strain: np.ndarray, total_mass: Number, distance: Number) -> np.ndarray:
    """Converts geometric strain into strain in meters.

    Args:
        strain (np.ndarray): Geometric strain.
        total_mass (Number): Total mass in solar masses.
        distance (Number): Distance in mega parsec.

    Returns:
        np.ndarray: Strain in meters.
    """
    return strain * total_mass * MTSUN_SI * C_SI / (distance * 1e6 * PC_SI)
