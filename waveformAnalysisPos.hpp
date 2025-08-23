#ifndef WAVEFORMANALYSISPOS_HPP
#define WAVEFORMANALYSISPOS_HPP

#include <string>
#include <vector>

struct Parameter {
  std::string label{};  // Parameter's name
  double min{};         // Min of range
  double max{};         // Max of range
};

// Choose area conversion factor
enum AreaConvFactor { Transm = 4200, Refl = 12500 };

// Pulse struct allows to define the properties of pulses in a waveform
struct Pulse {
  double startTime{};            // Overall start time of the pulse
  double endTime{};              // Overall end time of the pulse
  double peakTime{};             // Overall time at which peak takes place
  double peakValue{};            // Value of the pulse's peak in ADC
  double riseTime{};             // Time it takes to go from 10% to 90%
  double FWHMTime{};             // Time correspondent with the FWHM of the peak
  double areaFractionTime{};     // Fractional area time
  double area{};                 // Area of pulse
  double negFracArea{};          // Negative / overall area fraction
  double negFrac{};              // Negative / overall samples
  std::vector<double> values{};  // Values defining the pulse
  std::vector<double> times{};   // Times defining the pulse
};

// Waveform class allows to perform waveform analysis
class WaveformAnalysisPos {
 public:
  WaveformAnalysisPos(std::vector<double> const &s = {}, double ts = {},
                      double sp = {});

  // Getters for waveform member data

  std::vector<double> const &getSamples() const;
  double getTimeStamp() const;
  double getSamplePeriod() const;
  double getBaseline() const;
  double getThreshold() const;
  std::vector<Pulse> const &getPulses() const;

  // Analyse a waveform by extracting its noise and its pulses
  void analyseWaveform();

  // Make the average of first 50 samples of a waveform
  void baseline(int nInitialSamples = 50);

  // Compute RMS to find waveform's noise
  double RMS(int nInitialSamples = 50);

  // Detect pulses
  void findPulses(double threshold = 80., double tolerance = 50.,
                  int minWidth = 20, int maxWidth = 20, int minSep = 1);

  // Find area of 1 pulse
  Pulse integratePulse(int pulseStart, int pulseEnd);

 private:
  std::vector<double> fSamples{};  // Vector of samples generating the waveform
  double fTimeStamp{};             // Overall timestamp of the waveform
  double fSamplePeriod{};          // Sampling period
  double fBaseline{};              // Baseline of the waveform
  double fThreshold{};             // Threshold for pulse detection
  std::vector<Pulse> fPulses{};    // Vector of pulses composing the waveform
};

#endif