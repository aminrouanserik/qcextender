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

phenom = Waveform.from_model("IMRPhenomD", [(2, 2)], **kwargs)
seob = Waveform.from_model("SEOBNRv4", [(2, 2)], **kwargs)

sim = DimensionlessWaveform.from_sim("SXS:BBH:3977")
sim10sm, freq, freqadj, time = sim.to_Waveform(20, 100, 10, 0, 0)

print(phenom.match(seob))
print(sim10sm.match(phenom))
# print(phenom.metadata)
# print(seob.metadata)
# print(sim.metadata)
# print(sim10sm.metadata)

# plt.plot(
#     phenom.time,
#     np.gradient(-np.unwrap(np.angle(phenom[2, 2])), phenom.time),
#     label="phenom",
# )
# plt.plot(
#     seob.time, np.gradient(-np.unwrap(np.angle(seob[2, 2])), seob.time), label="seob"
# )
# plt.plot(time, freq, label="sim")
# plt.legend()
# plt.show()

plt.plot(phenom.time, phenom[2, 2])
plt.plot(seob.time, seob[2, 2])
plt.plot(sim10sm.time, sim10sm[2, 2])
plt.show()
