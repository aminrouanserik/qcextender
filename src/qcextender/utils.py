import numpy as np
import lal
import lalsimulation as ls


def lal_waves(
    approximant,
    mass1,
    mass2,
    spin1,
    spin2,
    distance,
    coa_phase,
    delta_t,
    f_lower,
    f_ref,
    mode,
):
    long_asc_nodes = 0.0
    eccentricity = 0.0
    mean_per_ano = 0.0
    lal_pars = None
    aprox = eval("ls." + str(approximant))

    generateTD = ls.SimInspiralTDModesFromPolarizations(
        mass1 * lal.MSUN_SI,
        mass2 * lal.MSUN_SI,
        spin1[0],
        spin1[1],
        spin1[2],
        spin2[0],
        spin2[1],
        spin2[2],
        distance * 1e6 * lal.PC_SI,
        coa_phase,
        long_asc_nodes,
        eccentricity,
        mean_per_ano,
        delta_t,
        f_lower,
        f_ref,
        lal_pars,
        aprox,
    )
    h22lal = ls.SphHarmTimeSeriesGetMode(generateTD, mode[0], mode[1])
    h22 = h22lal.data.data
    time = np.arange(h22lal.data.length, dtype=float) * h22lal.deltaT
    return time, h22
