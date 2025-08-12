# Params info
## Oscilloscope data for $30^{\circ}$
| Region     | Inc [PE/ns] | Transm [PE/ns] | Refl [PE/ns] | Prob_T | Prob_R |
|------------|-------------|----------------|--------------|--------|--------|
| PreTrig    | 1.897       | 1.700          | 2.387        | 0.896  | 1.258  |
| Trig       | 1055.665    | 229.932        | 380.482      | 0.218  | 0.360  |
| PostTrig1  | 472.744     | 117.550        | 150.555      | 0.249  | 0.318  |
| PostTrig2  | 254.758     | 84.519         | 74.242       | 0.332  | 0.291  |

### [4LayersA.txt](./data/4LayersA.txt)
```
  // Bin for hPhotoelectron = 60

  // Create fit function for summed pulses
  TF1 *fGaus = new TF1("fGaus", "gaus", 0.7, 1.5);

  // Make the average of first 1100 samples of a waveform
  void baseline(int nInitialSamples = 1100);

  // Compute RMS to find waveform's noise
  double RMS(int nInitialSamples = 1100);

  // Detect pulses
  void findPulses(double threshold = 19., double tolerance = 7.,
                  int minWidth = 20, int maxWidth = 20, int minSep = 10);

  // Area conversion factor (current assumption is 1 PE = 55 mV*ns)
  double const areaConvFactor{55.};

  // Below gaussian fit on "Pulse sum" graph is used
  double const triggerStart{1.1410014154 - 2 * 0.2461653408};
  double const triggerEnd{1.1410014154 + 2 * 0.2461653408};

  // Variable for analysis in pre-trigger region in 1 single file
  int pulseCounterPreTriggerRegion{};
  double totPreTrigArea{};
  double numPreTrigPE{};
  double const preTriggerStart{0.1};
  double const preTriggerEnd{0.7};

  // Variable for analysis in post-trigger region 1 in 1 single file
  int pulseCounterPostTriggerRegion1{};
  double totPostTrigArea1{};
  double numPostTrigPE1{};
  double const postTriggerStart1{1.7};
  double const postTriggerEnd1{2.3};

  // Variable for analysis in post-trigger region 2 in 1 single file
  int pulseCounterPostTriggerRegion2{};
  double totPostTrigArea2{};
  double numPostTrigPE2{};
  double const postTriggerStart2{3.};
  double const postTriggerEnd2{3.6};
```

### [incT.txt](./data/incT.txt)
```
  // Bin for hPhotoelectron = 60

  // Create fit function for summed pulses
  TF1 *fGaus = new TF1("fGaus", "gaus", 0.5, 0.8);
  
  // Creating TFile
  TFile *file1 = new TFile("./rootFiles/wA0.root", "RECREATE");

  // Make the average of first 170 samples of a waveform
  void baseline(int nInitialSamples = 170);

  // Compute RMS to find waveform's noise
  double RMS(int nInitialSamples = 170);

  // Detect pulses
  void findPulses(double threshold = 19., double tolerance = 7.,
                  int minWidth = 20, int maxWidth = 20, int minSep = 10);

  // Area conversion factor (current assumption is 1 PE = 37 mV*ns)
  double const areaConvFactor{37.};

  // Below gaussian fit on "Pulse sum" graph is used
  double const triggerStart{0.5742169474 - 2 * 0.0945427226};
  double const triggerEnd{0.5742169474 + 2 * 0.0945427226};

  // Variable for analysis in pre-trigger region in 1 single file
  int pulseCounterPreTriggerRegion{};
  double totPreTrigArea{};
  double numPreTrigPE{};
  double const preTriggerStart{0.};
  double const preTriggerEnd{0.45};

  // Variable for analysis in post-trigger region 1 in 1 single file
  int pulseCounterPostTriggerRegion1{};
  double totPostTrigArea1{};
  double numPostTrigPE1{};
  double const postTriggerStart1{0.8};
  double const postTriggerEnd1{1.2};

  // Variable for analysis in post-trigger region 2 in 1 single file
  int pulseCounterPostTriggerRegion2{};
  double totPostTrigArea2{};
  double numPostTrigPE2{};
  double const postTriggerStart2{1.4};
  double const postTriggerEnd2{1.8};
```

### [incR.txt](./data/incR.txt)
```
  // Bin for hPhotoelectron = 25

  // Create fit function for summed pulses
  TF1 *fGaus = new TF1("fGaus", "gaus", 0.5, 0.8);

  // Creating TFile
  TFile *file1 = new TFile("./rootFiles/wA1.root", "RECREATE");

  // Make the average of first 170 samples of a waveform
  void baseline(int nInitialSamples = 170);

  // Compute RMS to find waveform's noise
  double RMS(int nInitialSamples = 170);

  // Detect pulses
  void findPulses(double threshold = 19., double tolerance = 7.,
                  int minWidth = 20, int maxWidth = 20, int minSep = 10);

  // Area conversion factor (current assumption is 1 PE = 37 mV*ns)
  double const areaConvFactor{60.};

  // Below gaussian fit on "Pulse sum" graph is used
  double const triggerStart{0.5751793487 - 2 * 0.1677945632};
  double const triggerEnd{0.5751793487 + 2 * 0.1677945632};

  // Variable for analysis in pre-trigger region in 1 single file
  int pulseCounterPreTriggerRegion{};
  double totPreTrigArea{};
  double numPreTrigPE{};
  double const preTriggerStart{0.};
  double const preTriggerEnd{0.5};

  // Variable for analysis in post-trigger region 1 in 1 single file
  int pulseCounterPostTriggerRegion1{};
  double totPostTrigArea1{};
  double numPostTrigPE1{};
  double const postTriggerStart1{0.85};
  double const postTriggerEnd1{1.25};

  // Variable for analysis in post-trigger region 2 in 1 single file
  int pulseCounterPostTriggerRegion2{};
  double totPostTrigArea2{};
  double numPostTrigPE2{};
  double const postTriggerStart2{1.45};
  double const postTriggerEnd2{1.85};
```

### [transm.txt](./data/transm.txt)
```
  // Bin for hPhotoelectron = 60

  // Create fit function for summed pulses
  TF1 *fGaus = new TF1("fGaus", "gaus", 0.5, 0.8);

  // Creating TFile
  TFile *file1 = new TFile("./rootFiles/wA2.root", "RECREATE");

  // Make the average of first 170 samples of a waveform
  void baseline(int nInitialSamples = 170);

  // Compute RMS to find waveform's noise
  double RMS(int nInitialSamples = 170);

  // Detect pulses
  void findPulses(double threshold = 19., double tolerance = 7.,
                  int minWidth = 20, int maxWidth = 20, int minSep = 10);

  // Area conversion factor (current assumption is 1 PE = 37 mV*ns)
  double const areaConvFactor{37.};

  // Below gaussian fit on "Pulse sum" graph is used
  double const triggerStart{0.5880322313 - 2 * 0.0927908366};
  double const triggerEnd{0.5880322313 + 2 * 0.0927908366};

  // Variable for analysis in pre-trigger region in 1 single file
  int pulseCounterPreTriggerRegion{};
  double totPreTrigArea{};
  double numPreTrigPE{};
  double const preTriggerStart{0.};
  double const preTriggerEnd{0.5};

  // Variable for analysis in post-trigger region 1 in 1 single file
  int pulseCounterPostTriggerRegion1{};
  double totPostTrigArea1{};
  double numPostTrigPE1{};
  double const postTriggerStart1{0.85};
  double const postTriggerEnd1{1.25};

  // Variable for analysis in post-trigger region 2 in 1 single file
  int pulseCounterPostTriggerRegion2{};
  double totPostTrigArea2{};
  double numPostTrigPE2{};
  double const postTriggerStart2{1.45};
  double const postTriggerEnd2{1.85};
```

### [refl.txt](./data/refl.txt)
```
  // Bin for hPhotoelectron = 25

  // Create fit function for summed pulses
  TF1 *fGaus = new TF1("fGaus", "gaus", 0.5, 0.8);

  // Creating TFile
  TFile *file1 = new TFile("./rootFiles/wA3.root", "RECREATE");

  // Make the average of first 170 samples of a waveform
  void baseline(int nInitialSamples = 170);

  // Compute RMS to find waveform's noise
  double RMS(int nInitialSamples = 170);

  // Detect pulses
  void findPulses(double threshold = 19., double tolerance = 7.,
                  int minWidth = 20, int maxWidth = 20, int minSep = 10);

  // Area conversion factor (current assumption is 1 PE = 37 mV*ns)
  double const areaConvFactor{60.};

  // Below gaussian fit on "Pulse sum" graph is used
  double const triggerStart{0.5745826154 - 2 * 0.1743512132};
  double const triggerEnd{0.5745826154 + 2 * 0.1743512132};

  // Variable for analysis in pre-trigger region in 1 single file
  int pulseCounterPreTriggerRegion{};
  double totPreTrigArea{};
  double numPreTrigPE{};
  double const preTriggerStart{0.};
  double const preTriggerEnd{0.5};

  // Variable for analysis in post-trigger region 1 in 1 single file
  int pulseCounterPostTriggerRegion1{};
  double totPostTrigArea1{};
  double numPostTrigPE1{};
  double const postTriggerStart1{0.95};
  double const postTriggerEnd1{1.35};

  // Variable for analysis in post-trigger region 2 in 1 single file
  int pulseCounterPostTriggerRegion2{};
  double totPostTrigArea2{};
  double numPostTrigPE2{};
  double const postTriggerStart2{1.45};
  double const postTriggerEnd2{1.85};
```