//
// Created by mreichmann on 04.09.17.
//

#ifndef TRACKINGTELESCOPE_UTILS_H
#define TRACKINGTELESCOPE_UTILS_H

#include <string>
#include <vector>
#include <sys/ioctl.h>
#include <iostream>

/** terminal output colors */
#define ERROR "\033[91m"
#define ENDC "\033[0m"

struct winsize;

namespace tel{

  std::string trim(const std::string &, std::string="\t\n\r\v ");
  std::vector<std::string> split(const std::string &, const char & = '\t');
  double distance(std::pair<float, float>, std::pair<float, float>);
  void critical(const std::string & msg);

  class ProgressBar {
  private:
    uint32_t nEvents;
    uint32_t currentEvent;
    bool useETA;
    struct winsize w;
    uint8_t barLength;
    uint16_t updateFrequency;
    clock_t lastTime;
    uint8_t nCycles;
    float timePerCycle;

  public:
    explicit ProgressBar(uint32_t, bool use_ETA=true, uint16_t update=100);
    ~ProgressBar() = default;
    void update(uint32_t=0);
    void averageTime();
    float getTime();
    void reset() { currentEvent = 0u; }
  };
}

#endif //TRACKINGTELESCOPE_UTILS_H