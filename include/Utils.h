//
// Created by mreichmann on 04.09.17.
//

#ifndef TRACKINGTELESCOPE_UTILS_H
#define TRACKINGTELESCOPE_UTILS_H

#include <string>
#include <vector>
#include <sys/ioctl.h>
#include <iostream>
#include <sstream>

/** terminal output colors */
#define ERROR "\033[91m"
#define ENDC "\033[0m"

struct winsize;

namespace tel{

  std::string trim(const std::string &, const std::string & chars="\t\n\r\v ");
  std::vector<std::string> split(const std::string &, const char & = '\t');
  double distance(std::pair<float, float>, std::pair<float, float>);
  void critical(const std::string & msg);
  void print_banner(const std::string& message, char seperator = '=', uint16_t max_lenght = 100);


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
    void setNEvents(uint32_t n_events) { nEvents = n_events; }
    ProgressBar & operator++();
  };

  template <typename Q>
  inline std::string to_string(const std::vector<Q> & vec){
    std::ostringstream s;
    s << "[";
    for (int word: vec){
      s << word;
      s << (word == vec.back() ? "" : ", ");
    }
    s << "]";
    return s.str();
  }
}

#endif //TRACKINGTELESCOPE_UTILS_H