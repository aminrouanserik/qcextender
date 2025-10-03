import matplotlib.pyplot as plt
from qcextender.waveform import Waveform
from qcextender.dimensionlesswaveform import DimensionlessWaveform

import numpy as np

kwargs = {
    "mass1": 50,
    "mass2": 50,
    "inclination": 0,
    "coa_phase": 0,
    "delta_t": 1.0 / 4096,
    "f_lower": 20,
    "f_ref": 25,  # Change to be specific to waveform model used, want to do that in the generation.
    "distance": 10,
}

kwargsseob = {
    "mass1": 50,
    "mass2": 50,
    "inclination": 0,
    "coa_phase": np.pi / 2,
    "delta_t": 1.0 / 4096,
    "f_lower": 20,
    "f_ref": 25,  # Change to be specific to waveform model used, want to do that in the generation.
    "distance": 10,
}

phenom = Waveform.from_model("IMRPhenomD", [(2, 2)], **kwargs)
seob = Waveform.from_model("SEOBNRv4", [(2, 2)], **kwargsseob)

sim = DimensionlessWaveform.from_sim("SXS:BBH:3977")
sim10sm = sim.to_Waveform(20, 100, 10, 0, 0)

print(phenom.match(seob))
print(sim10sm.match(phenom))
# print(phenom.metadata)
# print(seob.metadata)
# print(sim.metadata)
print(sim10sm.metadata)

phenomfreq = phenom.freq()
seobfreq = seob.freq()
sim10smfreq = sim10sm.freq()

plt.plot(phenomfreq.sample_frequencies, phenomfreq.real(), label="IMRPhenomD")
plt.plot(seobfreq.sample_frequencies, seobfreq.real(), label="SEOBNRv4")
plt.plot(sim10smfreq.sample_frequencies, sim10smfreq.real(), label="SXS:BBH:3977")
plt.ylabel("Strain (m)")
plt.xlabel("Frequency (Hz)")
plt.legend()
plt.show()

plt.plot(phenom.time, phenom.phase(), label="IMRPhenomD")
plt.plot(seob.time, seob.phase(), label="SEOBNRv4")
plt.plot(sim10sm.time, sim10sm.phase(), label="SXS:BBH:3977")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude(m)")
plt.legend()
plt.show()

plt.plot(phenom.time, phenom.amp(), label="IMRPhenomD")
plt.plot(seob.time, seob.amp(), label="SEOBNRv4")
plt.plot(sim10sm.time, sim10sm.amp(), label="SXS:BBH:3977")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude(m)")
plt.legend()
plt.show()

plt.plot(phenom.time, phenom[2, 2], label="IMRPhenomD")
plt.plot(seob.time, seob[2, 2], label="SEOBNRv4")
plt.plot(sim10sm.time, sim10sm[2, 2], label="SXS:BBH:3977")
plt.xlabel("Time (s)")
plt.ylabel("Strain (m)")
plt.legend()
plt.show()
