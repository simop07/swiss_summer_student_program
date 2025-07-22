#include "waveformAnalysis.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>

WaveformAnalysis::WaveformAnalysis(std::vector<double> const &s, double ts,
                                   double sp)
    : fSamples{s}, fTimestamp{ts}, fSamplePeriod{sp} {
  WaveformAnalysis::analyse();  // When a waveform is found, analyse it
}

std::vector<double> const &WaveformAnalysis::getSamples() const {
  return fSamples;
}
double WaveformAnalysis::getTimestamp() const { return fTimestamp; }
double WaveformAnalysis::getSamplePeriod() const { return fSamplePeriod; }
double WaveformAnalysis::getBaseline() const { return fBaseline; }
std::vector<Pulse> const &WaveformAnalysis::getPulses() const {
  return fPulses;
}

void WaveformAnalysis::analyse() {
  WaveformAnalysis::computeBaseline();
  WaveformAnalysis::findAndIntegratePulse();
}

void WaveformAnalysis::computeBaseline(int nBaselineSamples) {
  fBaseline = std::accumulate(fSamples.begin(),
                              fSamples.begin() + nBaselineSamples, 0.) /
              nBaselineSamples;
}

double WaveformAnalysis::computeRMS(int nBaselineSamples) {
  double mean = fBaseline;
  double sum{};
  for (int i = 0; i < nBaselineSamples; i++) {
    double delta = fSamples[i] - mean;
    sum += delta * delta;  // Square the differences
  }
  return std::sqrt(sum / nBaselineSamples);
}

void WaveformAnalysis::findAndIntegratePulse(double tresholdFactor,
                                             double tolerance, int minWidth,
                                             int maxWidth, int minSeparation) {
  double const noiseRMS = computeRMS(50);

  // Threshold for peak detection
  double const pulseThreshold = fBaseline + tresholdFactor * noiseRMS;

  // Used to find the lower value of the interval for pulseStart and pulseEnd
  double const lowerLimit = fBaseline - tolerance * noiseRMS;

  // Used to find the higher value of the interval for pulseStart and pulseEnd
  double const upperLimit = fBaseline + tolerance * noiseRMS;

  bool inPulse = false;  // Am I in pulse?
  int pulseStart = 0;
  int pulseEnd = 0;

  // Min separation between pulseEnd of previous pulse and next pulse peak
  int lastPulseEnd = -minSeparation;

  for (int i = 0; i < static_cast<int>(fSamples.size()) - 2; i++) {
    // Detection of maximum value
    if (!inPulse && i - lastPulseEnd >= minSeparation &&
        fSamples[i] > pulseThreshold && fSamples[i] > fSamples[i - 1] &&
        fSamples[i] > fSamples[i + 1]) {
      inPulse = true;

      // Try to find startPulse near baseline before the peak
      bool foundStart = false;  // Have I found startPulse?
      for (int j = i; j >= std::max(0, i - minWidth); --j) {
        if (fSamples[j] >= lowerLimit && fSamples[j] <= upperLimit) {
          pulseStart = j;
          foundStart = true;
          break;
        }
      }
    }

    if (inPulse) {
      // Try to find endPulse near baseline after the peak
      bool foundEnd = false;  // Have I found endPulse?
      for (int j = i;
           j <= std::min(static_cast<int>(fSamples.size()) - 1, i + maxWidth);
           ++j) {
        if (fSamples[j] >= lowerLimit && fSamples[j] <= upperLimit) {
          pulseEnd = j;
          foundEnd = true;
          break;
        }
      }

      // Save pulse
      fPulses.push_back(integratePulse(pulseStart, pulseEnd));
      lastPulseEnd = pulseEnd;
      inPulse = false;
    }
  }

  // Limit case: still in pulse at the end
  if (inPulse) {
    pulseEnd = fSamples.size() - 1;
    fPulses.push_back(integratePulse(pulseStart, pulseEnd));
  }

  // Print waveform's properties
  std::cout << std::fixed << std::setprecision(1);  // Use 1 decimal digit
  std::cout << "\n\n=== PROPERTIES OF THE FOLLOWING WF ===\n";
  std::cout << "RMS           = " << noiseRMS << '\n';
  std::cout << "pulseThreshold= " << pulseThreshold << '\n';
  std::cout << "lowerLimit    = " << lowerLimit << '\n';
  std::cout << "upperLimit    = " << upperLimit << '\n';
}

// Find area of 1 pulse
Pulse WaveformAnalysis::integratePulse(int pulseStart, int pulseEnd) {
  int maxPulseId = pulseStart;
  double maxPulse = fSamples[pulseStart];
  double area{};

  // Find max of pulse
  for (int i = pulseStart; i <= pulseEnd; i++) {
    auto val = fSamples[i];
    if (val > maxPulse) {
      maxPulse = val;
      maxPulseId = i;
    }
    area += fSamples[i] - fBaseline;  // Remove baseline
  }
  // Right dimensions of area in ADC*ns
  area *= fSamplePeriod;

  return Pulse{fTimestamp + pulseStart * fSamplePeriod,
               fTimestamp + pulseEnd * fSamplePeriod,
               fTimestamp + maxPulseId * fSamplePeriod, maxPulse, area};
}