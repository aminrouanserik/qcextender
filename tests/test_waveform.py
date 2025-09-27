import matplotlib.pyplot as plt
from qcextender.waveform import Waveform

kwargs = {
    "mass1": 10,
    "mass2": 10,
    "inclination": 0,
    "coa_phase": 0,
    "delta_t": 1.0 / 4096,
    "f_lower": 20,
    "distance": 10,
}

phenom = Waveform.from_model("IMRPhenomD", [(2, 2), (3, 3)], **kwargs)
seob = Waveform.from_model("SEOBNRv4", [(2, 2)], **kwargs)

sim = Waveform.from_sim("SXS:BBH:3977", [(2, 2)])

print(phenom.metadata)
print(seob.metadata)
print(sim.metadata)

plt.plot(phenom.time, phenom[3, 3])
plt.plot(seob.time, seob[2, 2])
plt.plot(sim.time, sim[2, 2])
plt.show()
