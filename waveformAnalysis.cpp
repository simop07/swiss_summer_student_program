#include "waveformAnalysis.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>

WaveformAnalysis::WaveformAnalysis(std::vector<double> const &s, double ts,
                                   double sp)
    : fSamples{s}, fTimeStamp{ts}, fSamplePeriod{sp} {
  WaveformAnalysis::analyseWaveform();  // When a waveform is found, analyse it
}

std::vector<double> const &WaveformAnalysis::getSamples() const {
  return fSamples;
}

double WaveformAnalysis::getTimeStamp() const { return fTimeStamp; }

double WaveformAnalysis::getSamplePeriod() const { return fSamplePeriod; }

double WaveformAnalysis::getBaseline() const { return fBaseline; }

std::vector<Pulse> const &WaveformAnalysis::getPulses() const {
  return fPulses;
}

void WaveformAnalysis::analyseWaveform() {
  WaveformAnalysis::baseline();
  WaveformAnalysis::findPulses();
}

void WaveformAnalysis::baseline(int nInitialSamples) {
  // Sum the first nInitialSamples values and make their average
  auto sumSamples =
      std::accumulate(fSamples.begin(), fSamples.begin() + nInitialSamples, 0.);

  // Save value to member data
  fBaseline = sumSamples / nInitialSamples;
}

double WaveformAnalysis::RMS(int nInitialSamples) {
  double averageBaseline{fBaseline};
  double sumRMS{};
  for (int i{}; i < nInitialSamples; ++i) {
    // Compute differences with the average baseline
    double diff = fSamples[i] - averageBaseline;

    // Compute squares of differences
    sumRMS += diff * diff;
  }
  return std::sqrt(sumRMS / nInitialSamples);
}

void WaveformAnalysis::findPulses(double threshold, double tolerance,
                                  int minWidth, int maxWidth, int minSep) {
  // Compute RMS on waveforms' baseline
  double const rms = RMS(50);

  // Threshold for peak detection
  double const pulseThreshold = fBaseline + threshold * rms;

  // Used to find the lower value of the interval for pulseStart and pulseEnd
  double const lowLimit = fBaseline - tolerance * rms;

  // Used to find the higher value of the interval for pulseStart and pulseEnd
  double const upLimit = fBaseline + tolerance * rms;

  bool inPulse{false};  // Am I in pulse?
  int pulseStart{};
  int pulseEnd{};

  // Min separation between pulseEnd of previous pulse and next pulse peak
  int prevPulseEnd = -minSep;

  for (int i{}; i < static_cast<int>(fSamples.size()) - 2; ++i) {
    // Detection of pulse's maximum value
    if (!inPulse && fSamples[i] > pulseThreshold &&
        i - prevPulseEnd >= minSep && fSamples[i] > fSamples[i - 1] &&
        fSamples[i] > fSamples[i + 1]) {
      inPulse = true;

      // Try to find startPulse near baseline before the peak
      bool foundStart{false};  // Have I found startPulse?
      for (int j{i}; j >= std::max(0, i - minWidth); --j) {
        if (fSamples[j] > lowLimit && fSamples[j] < upLimit) {
          foundStart = true;
          pulseStart = j;
          break;
        }
      }
    }

    if (inPulse) {
      // Try to find endPulse near baseline after the peak
      bool foundEnd{false};  // Have I found endPulse?
      for (int j{i};
           j <= std::min(static_cast<int>(fSamples.size()) - 1, i + maxWidth);
           ++j) {
        if (fSamples[j] > lowLimit && fSamples[j] < upLimit) {
          foundEnd = true;
          pulseEnd = j;
          break;
        }
      }

      // Build pulse and save it into a vector
      fPulses.push_back(integratePulse(pulseStart, pulseEnd));
      prevPulseEnd = pulseEnd;
      inPulse = false;
    }
  }

  // Limit case where we are still in pulse at the end of the waveform
  if (inPulse) {
    // The endPulse becomes the last element of samples vector
    pulseEnd = static_cast<int>(fSamples.size()) - 1;

    // Build this pulse whatsoever
    fPulses.push_back(integratePulse(pulseStart, pulseEnd));
  }

  // Print waveform's properties
  std::cout << std::fixed << std::setprecision(2);  // Use 2 decimal digit
  std::cout << "\n\n*** PROPERTIES OF THE FOLLOWING WF ***\n";
  std::cout << "RMS value for noise             = " << rms << '\n';
  std::cout << "Threshold for pulse detection   = " << pulseThreshold << '\n';
  std::cout << "Lower limit for pulse endpoints = " << lowLimit << '\n';
  std::cout << "Upper limit for pulse endpoints = " << upLimit << '\n';
}

// Find area of 1 pulse
Pulse WaveformAnalysis::integratePulse(int pulseStart, int pulseEnd) {
  int maxPulseIndex = pulseStart;          // Index where maximum is verified
  double maxPulse = fSamples[pulseStart];  // Maximum value of the pulse
  double pulseArea{};

  // Find maximum value of the pulse
  for (int i = pulseStart; i <= pulseEnd; ++i) {
    auto currentVal = fSamples[i];
    if (currentVal > maxPulse) {
      maxPulse = currentVal;
      maxPulseIndex = i;
    }
    // Remove baseline to get correct area
    pulseArea += fSamples[i] - fBaseline;
  }

  // Right dimensions of area are in ADC*ns -> multiply by samplePeriod
  pulseArea *= fSamplePeriod;

  // Create pulse object
  return Pulse{fTimeStamp + pulseStart * fSamplePeriod,
               fTimeStamp + pulseEnd * fSamplePeriod,
               fTimeStamp + maxPulseIndex * fSamplePeriod, maxPulse, pulseArea};
}