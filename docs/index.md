# qcextender

**qcextender** is a Python package for extending and unifying gravitational waveform generation across libraries such as PyCBC, lalsimulation, and SXS.  
It provides a shared interface for handling metadata, performing unit conversions between dimensionless and physical (SI) representations, and managing consistent waveform transformations and match calculations.

The goal is to make it straightforward to load, compare, and analyze waveforms from different sources under a consistent, physically meaningful framework.

---

## Installation

### Using [uv](https://docs.astral.sh/uv/)

```bash
uv pip install git+https://github.com/aminrouanserik/qcextender.git
```

## Quick Example

```python
import matplotlib.pyplot as plt
from qcextender.waveform import Waveform
from qcextender.dimensionlesswaveform import DimensionlessWaveform

# Waveform parameters
mass1 = mass2 = 25
distance = 10
f_lower = 20
inclination, coa_phase = 0, 0
kwargs = {
    "mass1": mass1,
    "mass2": mass2,
    "inclination": inclination,
    "coa_phase": coa_phase,
    "delta_t": 1.0 / 4096,
    "f_lower": f_lower,
    "f_ref": 25,
    "distance": distance,
}

phenom = Waveform.from_model("IMRPhenomD", [(2, 2)], **kwargs)

# Load simulation and scale to SI units
sim = DimensionlessWaveform.from_sim("SXS:BBH:1155")
sim10sm = sim.to_Waveform(f_lower, mass1 + mass2, distance, inclination, coa_phase)

# Plot strain
plt.plot(phenom.time, phenom[2, 2], label="IMRPhenomD")
plt.plot(sim10sm.time, sim10sm[2, 2], label="SXS:BBH:1155")
plt.xlabel("Time (s)")
plt.ylabel("Strain (m)")
plt.legend()
plt.show()

```