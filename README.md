## qcextender

A Python package for extending and unifying gravitational waveform generation 
(PyCBC, SEOBNR, SXS, etc.) with consistent metadata, conversions between dimensionless and 
dimensioned representations, and standardized match calculations.

### Development Roadmap
Below is the current development status (as of v0.1.1):
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
    - [x] Shift amplitude ***and*** frequency (Currently shifting/cutting based on the phase (not a bad idea but needs correction))
1. **Adjust match calculation**
    - [x] Correct SXS delta_t
    - [x] Resize to single delta_t
    - [x] Calculate full waveform first
1. **Add eccentricity to Waveform**
    - [x] Returns waveform with changed eccentricity
1. **Clean up**
    *Refactor*
        - [x] Move (single-line) helper functions to utils (or units.py)
    *Docs*
        - [x] Add comments where necessary
        - [x] Adjust docstrings, add longer explanations next to single-line summaries
        - [ ] mkdocs full documentation

**Known Issues**
    - [ ] Change model and manually decompose into modes
    - [ ] Minor Metadata rework 
    - [ ] Clean SXS waveforms
    - [ ] Fix f_ref of other waveforms to SXS metadata
    - [ ] Check required keys, keyerror is vague
    - [ ] Tailored match caluclations