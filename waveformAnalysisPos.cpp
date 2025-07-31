#include "waveformAnalysisPos.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>

WaveformAnalysisPos::WaveformAnalysisPos(std::vector<double> const &s,
                                         double ts, double sp)
    : fSamples{s}, fTimeStamp{ts}, fSamplePeriod{sp} {
  WaveformAnalysisPos::analyseWaveform();  // When a waveform is found, analyse
                                           // it
}

std::vector<double> const &WaveformAnalysisPos::getSamples() const {
  return fSamples;
}

double WaveformAnalysisPos::getTimeStamp() const { return fTimeStamp; }

double WaveformAnalysisPos::getSamplePeriod() const { return fSamplePeriod; }

double WaveformAnalysisPos::getBaseline() const { return fBaseline; }

std::vector<Pulse> const &WaveformAnalysisPos::getPulses() const {
  return fPulses;
}

void WaveformAnalysisPos::analyseWaveform() {
  WaveformAnalysisPos::baseline();
  WaveformAnalysisPos::findPulses();
}

void WaveformAnalysisPos::baseline(int nInitialSamples) {
  // Sum the first nInitialSamples values and make their average
  auto sumSamples =
      std::accumulate(fSamples.begin(), fSamples.begin() + nInitialSamples, 0.);

  // Save value to member data
  fBaseline = sumSamples / nInitialSamples;
}

double WaveformAnalysisPos::RMS(int nInitialSamples) {
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

void WaveformAnalysisPos::findPulses(double threshold, double tolerance,
                                     int minWidth, int maxWidth, int minSep) {
  // Compute RMS on waveforms' baseline
  double const rms = RMS(50);

  // Threshold for peak detection
  double const pulseThreshold = fBaseline + threshold * rms;

  // Used to find the lower value of the interval for pulseStart and pulseEnd
  double const lowLimit = fBaseline - tolerance * rms;

  // Used to find the higher value of the interval for pulseStart and pulseEnd
  double const upLimit = fBaseline + tolerance * rms;

  int pulseStart{};
  int pulseEnd{};

  // Min separation between pulseEnd of previous pulse and next pulse peak
  int prevPulseEnd = -minSep;

  for (int i{1}; i < static_cast<int>(fSamples.size()) - 2; ++i) {
    // Detection of pulse's maximum value
    if (fSamples[i] > pulseThreshold && i - prevPulseEnd >= minSep &&
        fSamples[i] > fSamples[i - 1] && fSamples[i] > fSamples[i + 1]) {
      bool foundStart{false};  // Have I found startPulse?
      bool foundEnd{false};    // Have I found endPulse?

      // Try to find startPulse near baseline before the peak
      for (int j{i - 1}; j >= std::max(0, i - minWidth); --j) {
        if (fSamples[j] > lowLimit && fSamples[j] < upLimit) {
          foundStart = true;
          pulseStart = j;
          break;
        }
      }

      // Try to find endPulse near baseline after the peak if Start is found
      if (foundStart) {
        for (int j{i + 1};
             j <= std::min(static_cast<int>(fSamples.size()) - 1, i + maxWidth);
             ++j) {
          if (fSamples[j] > lowLimit && fSamples[j] < upLimit) {
            foundEnd = true;
            pulseEnd = j;
            break;
          } else if (j == static_cast<int>(fSamples.size()) - 1) {
            foundEnd = true;
            pulseEnd = j;
          }
        }
      }

      // Build pulse and save it into a vector
      if (foundStart && foundEnd && pulseStart < pulseEnd) {
        fPulses.push_back(integratePulse(pulseStart, pulseEnd));
        prevPulseEnd = pulseEnd;
      }
    }
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
Pulse WaveformAnalysisPos::integratePulse(int pulseStart, int pulseEnd) {
  // Find pulse maximum
  auto itMax = std::max_element(fSamples.begin() + pulseStart,
                                std::next(fSamples.begin() + pulseEnd));

  // Find index corresponding to the maximum
  auto maxPulseIndex = std::distance(fSamples.begin(), itMax);

  // Define peak to peak value
  double peakToPeak{std::abs(*itMax - fBaseline)};

  // Compute rise time
  auto itRiseBegin = std::find_if(
      fSamples.begin() + pulseStart, std::next(fSamples.begin() + pulseEnd),
      [&](double sample) { return sample >= (0.1 * peakToPeak + fBaseline); });
  auto itRiseEnd = std::find_if(
      fSamples.begin() + pulseStart, std::next(fSamples.begin() + pulseEnd),
      [&](double sample) { return sample >= (0.9 * peakToPeak + fBaseline); });
  // Find indices corresponding to rise time
  auto itRiseBeginIndex = std::distance(fSamples.begin(), itRiseBegin);
  auto itRiseEndIndex = std::distance(fSamples.begin(), itRiseEnd);
  double riseTime =
      static_cast<double>(itRiseEndIndex - itRiseBeginIndex) * fSamplePeriod;

  // Compute FWHM
  auto itFWHMBegin = std::find_if(
      fSamples.begin() + pulseStart, std::next(itMax),
      [&](double sample) { return sample >= (0.5 * peakToPeak + fBaseline); });
  auto itFWHMEnd = std::find_if(
      itMax, std::next(fSamples.begin() + pulseEnd),
      [&](double sample) { return sample <= (0.5 * peakToPeak + fBaseline); });
  auto itFWHMBeginIndex = std::distance(fSamples.begin(), itFWHMBegin);
  auto itFWHMEndIndex = std::distance(fSamples.begin(), itFWHMEnd);
  double FWHMTime =
      static_cast<double>(itFWHMEndIndex - itFWHMBeginIndex) * fSamplePeriod;

  // Compute sum of samples within a pulse
  auto sumSamples = std::accumulate(fSamples.begin() + pulseStart,
                                    std::next(fSamples.begin() + pulseEnd), 0.);

  // Compute sample average within a pulse
  auto avg = sumSamples / (pulseEnd - pulseStart + 1);

  // Subtract baseline in order to compute effective area
  auto pulseArea = sumSamples - (pulseEnd - pulseStart + 1) * fBaseline;

  // Find pulse which is closest to average
  int closestToAvgIndex = pulseStart;
  double minimalDifference = std::abs(fSamples[pulseStart] - avg);
  for (int i = pulseStart; i <= pulseEnd; ++i) {
    double currentVal = fSamples[i];
    auto diff = std::abs(currentVal - avg);
    if (diff < minimalDifference) {
      closestToAvgIndex = i;
      minimalDifference = diff;
    }
  }

  // Compute fractional area time (time width at which 90% of pulse area is
  // found, starting at the average point of the pulse)
  double pulseAreaPartial{};
  int rightTimeIndex{closestToAvgIndex};
  int leftTimeIndex{closestToAvgIndex - 1};
  while (pulseAreaPartial <= 0.9 * pulseArea) {
    if (rightTimeIndex <= pulseEnd) {
      pulseAreaPartial += (fSamples[rightTimeIndex] - fBaseline);
      ++rightTimeIndex;
    }
    if (leftTimeIndex >= pulseStart) {
      pulseAreaPartial += (fSamples[leftTimeIndex] - fBaseline);
      --leftTimeIndex;
    }
  }
  double fractionalAreaTime =
      (rightTimeIndex - leftTimeIndex - 1) * fSamplePeriod;

  // Using correct conversion for area
  pulseArea *= fSamplePeriod;

  // Create pulse object
  return Pulse{fTimeStamp + pulseStart * fSamplePeriod,
               fTimeStamp + pulseEnd * fSamplePeriod,
               fTimeStamp + maxPulseIndex * fSamplePeriod,
               *itMax,
               riseTime,
               FWHMTime,
               fractionalAreaTime,
               std::abs(pulseArea)};
}