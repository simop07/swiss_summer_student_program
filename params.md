# Params info

## Digitised data for $45^{\circ}$
### $1$ PTFE layer
| Region     | Inc_T [PE/ns] | Inc_R [PE/ns] | Transm [PE/ns] | Refl [PE/ns] | Prob_T | Prob_R |
|------------|---------------|---------------|----------------|--------------|--------|--------|
| PreTrig    | 0.000         | 0.000         | 0.000          | 0.000        | -nan   | -nan   |
| Trig       | 124.047       | 28.114        | 26.600         | 21.339       | 0.214  | 0.172  |
| PostTrig1  | 28.011        | 7.465         | 7.117          | 6.520        | 0.254  | 0.233  |
| PostTrig2  | 24.283        | 5.056         | 6.338          | 4.176        | 0.261  | 0.172  |

### $2$ PTFE layer
| Region     | Inc_T [PE/ns] | Inc_R [PE/ns] | Transm [PE/ns] | Refl [PE/ns] | Prob_T | Prob_R |
|------------|---------------|---------------|----------------|--------------|--------|--------|
| PreTrig    | 0.000         | 0.000         | 0.000          | 0.000        | -nan   | -nan   |
| Trig       | 124.047       | 28.114        | 15.240         | 25.203       | 0.123  | 0.203  |
| PostTrig1  | 28.011        | 7.465         | 4.719          | 6.952        | 0.168  | 0.248  |
| PostTrig2  | 24.283        | 5.056         | 4.084          | 6.530        | 0.168  | 0.269  |

### $3$ PTFE layer
| Region     | Inc_T [PE/ns] | Inc_R [PE/ns] | Transm [PE/ns] | Refl [PE/ns] | Prob_T | Prob_R |
|------------|---------------|---------------|----------------|--------------|--------|--------|
| PreTrig    | 0.000         | 0.000         | 0.000          | 0.000        | -nan   | -nan   |
| Trig       | 124.047       | 28.114        | 7.968          | 26.945       | 0.064  | 0.217  |
| PostTrig1  | 28.011        | 7.465         | 2.470          | 3.353        | 0.088  | 0.120  |
| PostTrig2  | 24.283        | 5.056         | 1.745          | 3.749        | 0.072  | 0.154  |

### Comparison table
| Thickness [mm] | Prob_T | Prob_R |
|----------------|--------|--------|
| 1.55           | 0.214  | 0.172  |
| 3.10           | 0.123  | 0.203  |
| 5.14           | 0.064  | 0.217  |

## Oscilloscope data for [4LayersA.txt](./data/4LayersA.txt)
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

## Oscilloscope data for $30^{\circ}$
| Region     | Inc [PE/ns] | Transm [PE/ns] | Refl [PE/ns] | Prob_T | Prob_R |
|------------|-------------|----------------|--------------|--------|--------|
| PreTrig    | 1.897       | 1.700          | 2.387        | 0.896  | 1.258  |
| Trig       | 1055.665    | 229.932        | 380.482      | 0.218  | 0.360  |
| PostTrig1  | 472.744     | 117.550        | 150.555      | 0.249  | 0.318  |
| PostTrig2  | 254.758     | 84.519         | 74.242       | 0.332  | 0.291  |

### [incT.txt](./data/incT.txt)
```
  // Bin for hPhotoelectron = 60

  // Create fit function for summed pulses
  TF1 *fGaus = new TF1("fGaus", "gaus", 0.5, 0.8);

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

## Oscilloscope data for $60^{\circ}$ ~ WORK IN PROGRESS
| Region     | Inc [PE/ns] | Transm [PE/ns] | Refl [PE/ns] | Prob_T | Prob_R |
|------------|-------------|----------------|--------------|--------|--------|
| PreTrig    | 1.147       | 0.000          | 3.370        | 0.000  | 2.938  |
| Trig       | 243.475     | 215.317        | 1043.009     | 0.884  | 4.284  |
| PostTrig1  | 174.684     | 101.659        | 439.060      | 0.582  | 2.513  |
| PostTrig2  | 111.794     | 80.781         | 189.913      | 0.723  | 1.699  |

### [inc60T.txt](./data/inc60T.txt)
```
  // Info hPhotoelectron
  TH1F *hPhotoElectrons =
      new TH1F("hPE", "Pulse area distribution; Area [PE]; Normalized counts",
               40, 0, 4.5);

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
  double const areaConvFactor{27.};

  // Below gaussian fit on "Pulse sum" graph is used
  double const triggerStart{0.5670213385 - 2 * 0.0938516294};
  double const triggerEnd{0.5670213385 + 2 * 0.0938516294};

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

### [transm60.txt](./data/transm60.txt)
```
  // Info hPhotoelectron
  TH1F *hPhotoElectrons =
      new TH1F("hPE", "Pulse area distribution; Area [PE]; Normalized counts",
               40, 0, 4.5);

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
  double const areaConvFactor{27.};

  // Below gaussian fit on "Pulse sum" graph is used
  double const triggerStart{0.5612464817 - 2 * 0.0749055716};
  double const triggerEnd{0.5612464817 + 2 * 0.0749055716};

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

### [refl60.txt](./data/refl60.txt)
```
  // Info for hPhotoelectron
  TH1F *hPhotoElectrons =
      new TH1F("hPE", "Pulse area distribution; Area [PE]; Normalized counts",
               40, 0, 4.5);

  // Create fit function for summed pulses
  TF1 *fGaus = new TF1("fGaus", "gaus", 0.5, 0.8);
  
  // Creating TFile
  TFile *file1 = new TFile("./rootFiles/wA60_2.root", "RECREATE");

  // Make the average of first 170 samples of a waveform
  void baseline(int nInitialSamples = 170);

  // Compute RMS to find waveform's noise
  double RMS(int nInitialSamples = 170);

  // Detect pulses
  void findPulses(double threshold = 19., double tolerance = 7.,
                  int minWidth = 20, int maxWidth = 20, int minSep = 10);

  // Area conversion factor (current assumption is 1 PE = 37 mV*ns)
  double const areaConvFactor{70.};

  // Below gaussian fit on "Pulse sum" graph is used
  double const triggerStart{0.5758587068 - 2 * 0.0889747690};
  double const triggerEnd{0.5758587068 + 2 * 0.0889747690};

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