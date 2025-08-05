#ifndef WAVEFORMANALYSISNEG_HPP
#define WAVEFORMANALYSISNEG_HPP

#include <vector>

struct Parameter {
  std::string label{};  // Parameter's name
  double min{};         // Min of range
  double max{};         // Max of range
};

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
  double posFracArea{};          // Positive / overall area fraction
  double posFrac{};              // Positive / overall samples
  std::vector<double> values{};  // Values defining the pulse
  std::vector<double> times{};   // Times defining the pulse
};

// Waveform class allows to perform waveform analysis
class WaveformAnalysisNeg {
 public:
  WaveformAnalysisNeg(std::vector<double> const &s = {}, double ts = {},
                      double sp = {});

  // Getters for waveform member data

  std::vector<double> const &getSamples() const;
  double getTimeStamp() const;
  double getSamplePeriod() const;
  double getBaseline() const;
  std::vector<Pulse> const &getPulses() const;

  // Analyse a waveform by extracting its noise and its pulses
  void analyseWaveform();

  // Make the average of first 1100 samples of a waveform
  void baseline(int nInitialSamples = 1100);

  // Compute RMS to find waveform's noise
  double RMS(int nInitialSamples = 1100);

  // Detect pulses
  void findPulses(double threshold = 20., double tolerance = 3.,
                  int minWidth = 30, int maxWidth = 30, int minSep = 10);

  // Find area of 1 pulse
  Pulse integratePulse(int pulseStart, int pulseEnd);

 private:
  std::vector<double> fSamples{};  // Vector of samples generating the waveform
  double fTimeStamp{};             // Overall timestamp of the waveform
  double fSamplePeriod{};          // Sampling period
  double fBaseline{};              // Baseline of the waveform
  std::vector<Pulse> fPulses{};    // Vector of pulses composing the waveform
};

#endif