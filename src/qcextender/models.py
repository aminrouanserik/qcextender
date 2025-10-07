"""Waveform generation interface using LALSimulation.

Provides a minimal, extendable interface for generating single-mode
gravitational waveforms in the time domain.

All quantities are returned in SI units. Masses are provided in solar
masses, distances in megaparsecs, and spins as dimensionless vectors.

Constants:
    PC_SI (float): Parsec in meters.
    MSUN_SI (float): Solar mass in kilograms.

Example:
    >>> from qcextender.waveform_generators import lal_mode
    >>> t, h = lal_mode(
    ...     approximant="IMRPhenomD",
    ...     mass1=30, mass2=25,
    ...     spin1=[0, 0, 0.5], spin2=[0, 0, 0.3],
    ...     distance=400, coa_phase=0,
    ...     delta_t=1/4096, f_lower=20, f_ref=20, mode=(2, 2)
    ... )
    >>> h.shape
    (16221,)
"""

import numpy as np
import lalsimulation as ls

PC_SI: float = 3.085677581491367e16
MSUN_SI: float = 1.9884098706980507e30


def lal_mode(
    approximant: str,
    mass1: float,
    mass2: float,
    spin1: np.ndarray | list[float],
    spin2: np.ndarray | list[float],
    distance: float,
    coa_phase: float,
    delta_t: float,
    f_lower: float,
    f_ref: float,
    mode: tuple[int, int],
) -> tuple[np.ndarray, np.ndarray]:
    """Generate a single gravitational-wave mode using LALSimulation. Wraps
    `      lalsimulation.SimInspiralTDModesFromPolarizations`.

        Args:
            approximant (str): Name of the LALSimulation approximant (e.g. "IMRPhenomD").
            mass1 (float): Primary mass [solar masses].
            mass2 (float): Secondary mass [solar masses].
            spin1 (array-like): Dimensionless spin vector (S1x, S1y, S1z).
            spin2 (array-like): Dimensionless spin vector (S2x, S2y, S2z).
            distance (float): Luminosity distance to source [megaparsecs].
            coa_phase (float): Coalescence phase [radians].
            delta_t (float): Sampling interval [seconds].
            f_lower (float): Lower frequency cutoff [Hz].
            f_ref (float): Reference frequency [Hz].
            mode (tuple[int, int]): Mode (l, m) to extract.

        Returns:
            tuple[np.ndarray, np.ndarray]:
                - time (np.ndarray): Time samples [seconds].
                - h (np.ndarray): Complex strain for the specified mode.
    """
    long_asc_nodes = 0.0
    eccentricity = 0.0
    mean_per_ano = 0.0
    lal_pars = None

    generateTD = ls.SimInspiralTDModesFromPolarizations(  # ls.SimInspiralTD
        mass1 * MSUN_SI,
        mass2 * MSUN_SI,
        spin1[0],
        spin1[1],
        spin1[2],
        spin2[0],
        spin2[1],
        spin2[2],
        distance * 1e6 * PC_SI,
        coa_phase,
        long_asc_nodes,
        eccentricity,
        mean_per_ano,
        delta_t,
        f_lower,
        f_ref,
        lal_pars,
        eval("ls." + str(approximant)),
    )

    hlal = ls.SphHarmTimeSeriesGetMode(generateTD, mode[0], mode[1])
    h = hlal.data.data
    time = np.arange(hlal.data.length) * hlal.deltaT

    return time, h


def lal_polarizations(
    approximant: str,
    mass1: float,
    mass2: float,
    spin1: np.ndarray | list[float],
    spin2: np.ndarray | list[float],
    distance: float,
    inclination: float,
    coa_phase: float,
    delta_t: float,
    f_lower: float,
    f_ref: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate a single gravitational-wave mode using LALSimulation. Wraps
    `      lalsimulation.SimInspiralTDModesFromPolarizations`.

        Args:
            approximant (str): Name of the LALSimulation approximant (e.g. "IMRPhenomD").
            mass1 (float): Primary mass [solar masses].
            mass2 (float): Secondary mass [solar masses].
            spin1 (array-like): Dimensionless spin vector (S1x, S1y, S1z).
            spin2 (array-like): Dimensionless spin vector (S2x, S2y, S2z).
            distance (float): Luminosity distance to source [megaparsecs].
            coa_phase (float): Coalescence phase [radians].
            delta_t (float): Sampling interval [seconds].
            f_lower (float): Lower frequency cutoff [Hz].
            f_ref (float): Reference frequency [Hz].
            mode (tuple[int, int]): Mode (l, m) to extract.

        Returns:
            tuple[np.ndarray, np.ndarray]:
                - time (np.ndarray): Time samples [seconds].
                - h (np.ndarray): Complex strain for the specified mode.
    """
    long_asc_nodes = 0.0
    eccentricity = 0.0
    mean_per_ano = 0.0
    lal_pars = None

    hp, hc = ls.SimInspiralTD(
        mass1 * MSUN_SI,
        mass2 * MSUN_SI,
        spin1[0],
        spin1[1],
        spin1[2],
        spin2[0],
        spin2[1],
        spin2[2],
        distance * 1e6 * PC_SI,
        inclination,
        coa_phase,
        long_asc_nodes,
        eccentricity,
        mean_per_ano,
        delta_t,
        f_lower,
        f_ref,
        lal_pars,
        eval("ls." + str(approximant)),
    )

    time = np.arange(hp.data.length) * hp.deltaT
    return hp.data.data, hc.data.data, time
