### **To-Do**
1. **Seperate Metadata**
    - [x] New Metadata file
    - [x] Avoid printing None features
    - [x] Dimensionless vs. dimensioned split on creation
    - [x] Change print order in __repr__()
    - [x] Conversion
1. **Split Waveform objects**
    - [x] BaseWaveform containing all common features, time-series
    - [x] SXS waveform, dimensionless, time-series
    - [x] Waveform, dimensioned, time-series
1. **Convert to Waveform**
    - [x] Add a feature turning SXS Waveforms into Waveforms when specifying mass, distance, etc.
    - [ ] Shift amplitude ***and*** frequency
1. **Adjust match calculation**
    - [ ] Turn into complex match (real and imaginary should be identical -> check)
1. **Add eccentricity to Waveform**
    - [ ] Insert a function allowing the eccentricity to be changed for a Waveform object

**Frequency scaling**
    - [x] Rewrite to save h_lm instead of hp - i hc
    - [x] Decompose h_lm into amplitude and phase (to_waveform)
    - [x] Cut at f_lower (to_waveform)
    - [x] Reconstruct (to_waveform)
    - For some reason does not work as intended. Plot (time vs freq for all 3 and see whats going on)
    - [ ] Add f_ref to Metadata and rework
    - [ ] Fix f_ref to SXS metadata
    - [ ] Investigate why (3, 3) mode returns None for IMRPhenomD