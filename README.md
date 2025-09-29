### **To-Do**
1. **Seperate Metadata**
    - [x] New Metadata file
    - [x] Avoid printing None features
    - [x] Dimensionless vs. dimensioned split on creation
    - [x] Change print order in __repr__()
    - [ ] Conversion and lock
1. **Split Waveform objects**
    - [x] BaseWaveform containing all common features, time-series
    - [x] SXS waveform, dimensionless, time-series
    - [x] Waveform, dimensioned, time-series
1. **Convert to Waveform**
    - [ ] Add a feature turning SXS Waveforms into Waveforms when specifying mass, distance, etc.
    - [ ] Shift amplitude ***and*** phase
1. **Adjust match calculation**
    - [ ] Turn into complex match?
1. **Add eccentricity to Waveform**
    - [ ] Insert a function allowing the eccentricity to be changed for a Waveform object