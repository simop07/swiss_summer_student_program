#ifndef WAVEFORMANALYSIS_HPP
#define WAVEFORMANALYSIS_HPP

#include <vector>

// Define constants
inline constexpr int nMinAnalysedRows{0};     // minimum index of analysed rows
inline constexpr int nMaxAnalysedRows{9962};  // maximum rows

// Pulse struct allows to define the properties of pulses in a waveform
struct Pulse {
  double startTime{};  // Overall start time of the pulse
  double endTime{};    // Overall end time of the pulse
  double peakTime{};   // Overall time at which peak takes place
  double peakValue{};  // Value of the pulse's peak in ADC
  double area{};       // Area of pulse
};

// Waveform class allows to perform waveform analysis
class WaveformAnalysis {
 public:
  WaveformAnalysis(std::vector<double> const &s, double ts, double sp);

  // Getters for waveform member data

  std::vector<double> const &getSamples() const;
  double getTimestamp() const;
  double getSamplePeriod() const;
  double getBaseline() const;
  std::vector<Pulse> const &getPulses() const;

  // Analyse a waveform by extracting its noise and its pulses
  void analyse();

  // Make the average of first 50 samples of a waveform
  void computeBaseline(int nBaselineSamples = 50);

  // Compute RMS to find waveform's noise
  double computeRMS(int nBaselineSamples = 50);

  // Detect pulses
  void findAndIntegratePulse(double tresholdFactor = 30.0,
                             double tolerance = 20, int minWidth = 3,
                             int maxWidth = 5, int minSeparation = 10);

  // Find area of 1 pulse
  Pulse integratePulse(int pulseStart, int pulseEnd);

 private:
  std::vector<double> fSamples{};  // Vector of samples forming the waveform
  double fTimestamp{};             // Overall timestamp of the waveform
  double fSamplePeriod{};          // Sampling period
  double fBaseline{};              // Baseline of the waveform
  std::vector<Pulse> fPulses{};    // Vector of pulses forming the waveform
};

#endif