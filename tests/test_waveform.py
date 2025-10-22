# DimensionlessWaveform will NOT be tested here anymore
# Once test routine is thought out, might include pytest

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from qcextender.waveform import Waveform

# from qcextender.dimensionlesswaveform import DimensionlessWaveforms


# # Setup constants
# spin1 = (0, 0, 0)
# spin2 = (0, 0, 0)
# q = 1
# distance, inclination, coa_phase = 10, 0, 0
# f_lower = 20

# # Parameter grid
# eccentricities = np.linspace(1e-8, 0.2, 20)
# total_masses = np.linspace(10, 200, 20)

# # Preallocate mismatch results
# mismatches = np.zeros((len(total_masses), len(eccentricities)))

# # Compute mismatches
# for i, total_mass in enumerate(total_masses):
#     for j, e in enumerate(eccentricities):
#         mass1 = q * total_mass / (q + 1)
#         mass2 = total_mass / (q + 1)

#         kwargs = {
#             "mass1": mass1,
#             "mass2": mass2,
#             "spin1": spin1,
#             "spin2": spin2,
#             "eccentricity": e,
#             "inclination": inclination,
#             "distance": distance,
#             "coa_phase": coa_phase,
#             "delta_t": 1.0 / 4096,
#             "f_ref": 20,
#             "f_lower": f_lower,
#         }

#         ecctd = Waveform.from_model("EccentricTD", [(2, 2)], **kwargs)
#         kwargs.pop("eccentricity")
#         kwargs.update({"amplitude0": 0, "phase0": 4})
#         tt4 = Waveform.from_model("TaylorT4", [(2, 2)], **kwargs)

#         mismatch = 1 - ecctd.match(tt4, f_lower)
#         mismatches[j, i] = max(mismatch, 1e-8)  # avoid zeros for log scale

# # Plot mismatch heatmap (log scale)
# plt.figure(figsize=(8, 6))
# plt.imshow(
#     mismatches,
#     origin="lower",
#     aspect="auto",
#     extent=[
#         total_masses.min(),
#         total_masses.max(),
#         eccentricities.min(),
#         eccentricities.max(),
#     ],
#     cmap="viridis",
#     norm=colors.LogNorm(vmin=np.nanmin(mismatches), vmax=np.nanmax(mismatches)),
# )
# plt.colorbar(label="Mismatch (log scale)")
# plt.xlabel(r"Total mass $M/M_\odot$")
# plt.ylabel(r"Eccentricity $e$")
# plt.title("Mismatch Heatmap: TaylorT4 vs EccentricTD")
# plt.tight_layout()
# plt.show()


# Setup constants
spin1 = (0, 0, 0)
spin2 = (0, 0, 0)
distance, inclination, coa_phase = 10, 0, 0
f_lower = 20

# Parameter grid
mass_ratios = np.linspace(1, 10, 20)
total_masses = np.linspace(20, 200, 20)

# Preallocate mismatch results
mismatches = np.zeros((len(total_masses), len(mass_ratios)))

# Compute mismatches
for i, total_mass in enumerate(total_masses):
    for j, q in enumerate(mass_ratios):
        mass1 = q * total_mass / (q + 1)
        mass2 = total_mass / (q + 1)

        kwargs = {
            "mass1": mass1,
            "mass2": mass2,
            "spin1": spin1,
            "spin2": spin2,
            "inclination": inclination,
            "distance": distance,
            "coa_phase": coa_phase,
            "delta_t": 1.0 / 20000,
            "f_ref": 20,
            "f_lower": f_lower,
            "eccentricity": 1e-8,
        }

        phen_ecc = Waveform.from_model("TaylorT4", [(2, 2)], **kwargs)
        phen_circ = Waveform.from_model("EccentricTD", [(2, 2)], **kwargs)

        mismatch = 1 - phen_ecc.match(phen_circ, f_lower)
        mismatches[i, j] = max(mismatch, 1e-8)  # avoid zeros for log scale

# Plot mismatch heatmap (log scale)
plt.figure(figsize=(8, 6))
plt.imshow(
    mismatches,
    origin="lower",
    aspect="auto",
    extent=[
        mass_ratios.min(),
        mass_ratios.max(),
        total_masses.min(),
        total_masses.max(),
    ],
    cmap="viridis",
    norm=colors.LogNorm(vmin=np.nanmin(mismatches), vmax=np.nanmax(mismatches)),
)
plt.colorbar(label="Mismatch (log scale)")
plt.xlabel("Mass Ratio (q = m₁/m₂)")
plt.ylabel("Total Mass [M☉]")
plt.title("Mismatch Heatmap: IMRPhenomD vs SEOBNRv4")
plt.tight_layout()
plt.show()

# phenom = Waveform.from_model("IMRPhenomD", [(2, 2)], **kwargs)
# seob = Waveform.from_model("SEOBNRv4", [(2, 2)], **kwargs)

# sim = DimensionlessWaveform.from_sim("SXS:BBH:1155")
# sim10sm = sim.to_Waveform(20, 50, 10, 0, 0)

# print(phenom.match(seob))
# print(sim10sm.match(phenom))

# print(phenom.metadata)
# print(seob.metadata)
# print(sim.metadata)
# print(sim10sm.metadata)

# phenomfreq = phenom.freq()
# seobfreq = seob.freq()
# sim10smfreq = sim10sm.freq()

# plt.plot(phenomfreq.sample_frequencies, phenomfreq.real(), label="IMRPhenomD")
# plt.plot(seobfreq.sample_frequencies, seobfreq.real(), label="SEOBNRv4")
# plt.plot(sim10smfreq.sample_frequencies, sim10smfreq.real(), label="SXS:BBH:1155")
# plt.ylabel("Strain (m)")
# plt.xlabel("Frequency (Hz)")
# plt.legend()
# plt.show()

# plt.plot(phenom.time, phenom.phase(), label="IMRPhenomD")
# plt.plot(seob.time, seob.phase(), label="SEOBNRv4")
# plt.plot(sim10sm.time, sim10sm.phase(), label="SXS:BBH:1155")
# plt.xlabel("Time (s)")
# plt.ylabel("Phase")
# plt.legend()
# plt.show()

# plt.plot(phenom.time, phenom.amp(), label="IMRPhenomD")
# plt.plot(seob.time, seob.amp(), label="SEOBNRv4")
# plt.plot(sim10sm.time, sim10sm.amp(), label="SXS:BBH:1155")
# plt.xlabel("Time (s)")
# plt.ylabel("Amplitude (m)")
# plt.legend()
# plt.show()

# plt.plot(phenom.time, phenom[2, 2], label="IMRPhenomD")
# plt.plot(seob.time, seob[2, 2], label="SEOBNRv4")
# plt.plot(sim10sm.time, sim10sm[2, 2], label="SXS:BBH:1155")
# plt.xlabel("Time (s)")
# plt.ylabel("Strain (m)")
# plt.legend()
# plt.show()

# Turn into test cases
# print(phenom.time[np.argmax(phenom.amp())] == 0)
# print(seob.time[np.argmax(seob.amp())] == 0)
# print(sim10sm.time[np.argmax(sim10sm.amp())] == 0)

# This or close to 0
# print(np.isclose(phenom.phase()[np.argmax(phenom.amp())] % np.pi, np.pi))
# print(np.isclose(seob.phase()[np.argmax(seob.amp())] % np.pi, np.pi))
# print(np.isclose(sim10sm.phase()[np.argmax(sim10sm.amp())] % np.pi, np.pi))
