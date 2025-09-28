### **To-Do**
1. **Seperate Metadata**
    - New Metadata file
    - Avoid printing None features
    - Dimensionless vs. dimensioned split
    - Frequency vs time split
1. **Split Waveform objects**
    - BaseWaveform containing all common features
    - SXS waveform, dimensionless time-series
    - Waveform, dimensioned time-series
1. **Convert to Waveform**
    - Add a feature turning SXS Waveforms into Waveforms when specifying mass, distance, etc.
    - Shift amplitude ***and*** phase
1. **Adjust match calculation**
    - Turn into complex match?
1. **Add eccentricity to Waveform**
    - Insert a function allowing the eccentricity to be changed for a Waveform object