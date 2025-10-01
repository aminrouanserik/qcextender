## qcextender

A Python package for extending and unifying gravitational waveform generation 
(PyCBC, SEOBNR, SXS, etc.) with consistent metadata, conversions between dimensionless and 
dimensioned representations, and standardized match calculations.

### Development Roadmap
Below is the current development status (as of v0.1.4):
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
    - [ ] Lower frequency bound
1. **Adjust match calculation**
    - [x] Correct SXS delta_t
    - [x] Resize to single delta_t
    - [x] Calculate full waveform first
1. **Add eccentricity to Waveform**
    - [ ] Add a function allowing the eccentricity to be changed for a Waveform object, using a function hook
1. **Clean up**
    *Refactor*
        - [ ] Move (single-line) helper functions to utils (or units.py)
        - [ ] Adjust Metadata class (add lock, f_ref and check conversion from DimensionlessWaveform -> Waveform)
    *Docs*
        - [ ] Add comments where necessary
        - [ ] Adjust docstrings, add longer explanations next to single-line summaries
        - [ ] mkdocs full documentation

**Issues**
    - [ ] Fix f_ref of other waveforms to SXS metadata
    - [ ] Investigate why (3, 3) mode returns None for IMRPhenomD

**Open Questions**
    - Decide whether to calculate the complex match since the real does not exactly equal the imaginary match.
    - Identify what needs to change in Metadata for modes instead of polarizations