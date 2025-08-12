# Params info
## Oscilloscope data
### [4LayersA.txt](./data/4LayersA.txt)
```
  // Make the average of first 1100 samples of a waveform
  void baseline(int nInitialSamples = 1100);

  // Compute RMS to find waveform's noise
  double RMS(int nInitialSamples = 1100);

  // Detect pulses
  void findPulses(double threshold = 19., double tolerance = 7.,
                  int minWidth = 20, int maxWidth = 20, int minSep = 10);

  // Area conversion factor (current assumption is 1 PE = 55 mV*ns)
  double const areaConvFactor{55.}; 
```

### [incT.txt](./data/incT.txt)
```
  // Make the average of first 170 samples of a waveform
  void baseline(int nInitialSamples = 170);

  // Compute RMS to find waveform's noise
  double RMS(int nInitialSamples = 170);

  // Detect pulses
  void findPulses(double threshold = 19., double tolerance = 7.,
                  int minWidth = 20, int maxWidth = 20, int minSep = 10);

  // Area conversion factor (current assumption is 1 PE = 37 mV*ns)
  double const areaConvFactor{37.};
```