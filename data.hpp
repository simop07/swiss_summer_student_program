#ifndef DATA_HPP
#define DATA_HPP

#include <vector>

struct RegionData {
  double PECounter{};  // Number of PE
  double deltaT{};     // Time extension
  double PEPulses{};   // Number of PE per pulse
};

struct PhotonData {
  RegionData preTrigger{};
  RegionData inTrigger{};
  RegionData postTrigger1{};
  RegionData postTrigger2{};
};

struct Point {
  double x{};
  double y{};
};

#endif