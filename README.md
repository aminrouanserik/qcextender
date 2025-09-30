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
    - [x] Shift amplitude ***and*** frequency
1. **Adjust match calculation**
    - [x] Correct SXS delta_t
    - [x] Resize to single delta_t
    - [ ] Match is now of polarizations, or calculate full wave first? (Open question)
    - [ ] If full wave, check real == imag (Note, only for full waveform so need to know previous answer)
1. **Add eccentricity to Waveform**
    - [ ] Add a function allowing the eccentricity to be changed for a Waveform object, using a function hook
1. **Clean up**
    - [ ] Move (single-line) helper functions to utils (or units.py)
    - [ ] Adjust Metadata class (add lock, f_ref and check conversion. Also polarizations, when are spherical harmonics expected in Metadata?)
    - [ ] Add comments where necessary
    - [ ] Adjust docstrings, add longer explanations next to single-line summaries
    - [ ] mkdocs full documentation

**Frequency scaling**
    - [x] Rewrite to save h_lm instead of hp - i hc
    - [x] Decompose h_lm into amplitude and phase (to_waveform)
    - [x] Cut at f_lower (to_waveform)
    - [x] Reconstruct (to_waveform)
    - [ ] Fix f_ref of other waveforms to SXS metadata
    - [ ] Investigate why (3, 3) mode returns None for IMRPhenomD