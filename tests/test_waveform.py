# DimensionlessWaveform will NOT be tested here anymore
# Once test routine is thought out, might include pytest

import numpy as np
import matplotlib.pyplot as plt
from qcextender.waveform import Waveform

from qcextender.dimensionlesswaveform import DimensionlessWaveform

kwargs = {
    "mass1": 25,
    "mass2": 25,
    "inclination": 0,
    "coa_phase": 0,
    "delta_t": 1.0 / 16300,
    "f_lower": 20,
    "f_ref": 25,
    "distance": 10,
}

phenom = Waveform.from_model("IMRPhenomD", [(2, 2)], **kwargs)
seob = Waveform.from_model("SEOBNRv4", [(2, 2)], **kwargs)

sim = DimensionlessWaveform.from_sim("SXS:BBH:1155")
sim10sm = sim.to_Waveform(20, 50, 10, 0, 0)

# print(phenom.match(seob))
# print(sim10sm.match(phenom))

# print(phenom.metadata)
# print(seob.metadata)
# print(sim.metadata)
# print(sim10sm.metadata)

phenomfreq = phenom.freq()
seobfreq = seob.freq()
sim10smfreq = sim10sm.freq()

plt.plot(phenomfreq.sample_frequencies, phenomfreq.real(), label="IMRPhenomD")
plt.plot(seobfreq.sample_frequencies, seobfreq.real(), label="SEOBNRv4")
plt.plot(sim10smfreq.sample_frequencies, sim10smfreq.real(), label="SXS:BBH:1155")
plt.ylabel("Strain (m)")
plt.xlabel("Frequency (Hz)")
plt.legend()
plt.show()

plt.plot(phenom.time, phenom.phase(), label="IMRPhenomD")
plt.plot(seob.time, seob.phase(), label="SEOBNRv4")
plt.plot(sim10sm.time, sim10sm.phase(), label="SXS:BBH:1155")
plt.xlabel("Time (s)")
plt.ylabel("Phase")
plt.legend()
plt.show()

plt.plot(phenom.time, phenom.amp(), label="IMRPhenomD")
plt.plot(seob.time, seob.amp(), label="SEOBNRv4")
plt.plot(sim10sm.time, sim10sm.amp(), label="SXS:BBH:1155")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude (m)")
plt.legend()
plt.show()

plt.plot(phenom.time, phenom[2, 2], label="IMRPhenomD")
plt.plot(seob.time, seob[2, 2], label="SEOBNRv4")
plt.plot(sim10sm.time, sim10sm[2, 2], label="SXS:BBH:1155")
plt.xlabel("Time (s)")
plt.ylabel("Strain (m)")
plt.legend()
plt.show()

# Turn into test cases
print(phenom.time[np.argmax(phenom.amp())] == 0)
print(seob.time[np.argmax(seob.amp())] == 0)
print(sim10sm.time[np.argmax(sim10sm.amp())] == 0)

# This or close to 0
print(np.isclose(phenom.phase()[np.argmax(phenom.amp())] % np.pi, np.pi))
print(np.isclose(seob.phase()[np.argmax(seob.amp())] % np.pi, np.pi))
print(np.isclose(sim10sm.phase()[np.argmax(sim10sm.amp())] % np.pi, np.pi))
