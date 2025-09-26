import matplotlib.pyplot as plt
from qcextender.waveform import Waveform

kwargs = {
    "mass1": 10,
    "mass2": 10,
    "spin1z": 0.9,
    "spin2z": 0.4,
    "inclination": 1.23,
    "coa_phase": 2.45,
    "delta_t": 1.0 / 4096,
    "f_lower": 40,
}

phenom = Waveform.from_model("IMRPhenomD", [(2, 2), (3, 3)], **kwargs)
seob = Waveform.from_model("SEOBNRv4", [(2, 2)], **kwargs)

print(type(phenom))

plt.plot(phenom.time, phenom.singlemode(3, 3))
plt.plot(seob.time, seob.singlemode())
plt.show()
