//
// Created by mreichmann on 04.09.17.
//

#include "Utils.h"
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

namespace tel {

  vector<string> split(const string & s, const char & deliminator) {

    stringstream ss(s);
    string item;
    vector<string> elements;
    while (getline(ss, item, deliminator))
      elements.push_back(item);
    return elements;
  }

  string trim(const string &str, const string & trim_characters) {

    size_t b = str.find_first_not_of(trim_characters);
    size_t e = str.find_last_not_of(trim_characters);
    if (b == std::string::npos || e == std::string::npos)
      return "";
    return string(str, b, e - b + 1);
  }

  double distance(std::pair<float, float> p1, std::pair<float, float> p2){

    return sqrt(pow(p1.first - p2.first, 2) + pow(p1.second - p2.second, 2));
  }

  ProgressBar::ProgressBar(uint32_t n_events, bool use_ETA, uint16_t update): nEvents(n_events), currentEvent(0), useETA(use_ETA), w(), updateFrequency(update), lastTime(clock()),
                                                                              nCycles(0), timePerCycle(0) {
    ioctl(0, TIOCGWINSZ, &w);
    uint8_t diff = not use_ETA ? uint8_t(26) : uint8_t(37);
    barLength = uint8_t(w.ws_col - diff);
  }

   void critical(const std::string & msg) { std::cerr << std::endl << ERROR << "CRITICAL: " << ENDC << msg << std::endl; }

  void ProgressBar::update(uint32_t event) {

    if (currentEvent == 0u)
      cout << endl;

    currentEvent = (event == 0u ? currentEvent + 1 : event);

    if (currentEvent % updateFrequency == 0 && (currentEvent != 0u)){
      stringstream ss;
      ss << "\rEvent: "  << setw(7) << setfill('0') << fixed << currentEvent;
      auto current_bar_length = uint8_t(float(currentEvent) / nEvents * barLength);
      ss << " [" << string(current_bar_length, '=') << ">" << string(barLength - current_bar_length, ' ') << "] ";
      ss << setprecision(2) << setw(6) << setfill('0') << fixed << float(currentEvent) / nEvents * 100 << "%";
      if (useETA) {
        averageTime();
        float tot = float(nEvents - currentEvent) / updateFrequency * timePerCycle;
        ss << setprecision(0) << " ETA: " << setw(2) << setfill('0') << tot / 60 << ":" << setw(2) << setfill('0') << tot - int(tot) / 60 * 60;
      }
      cout << ss.str() << flush;
    }
    if (currentEvent == nEvents) {
      cout << "\rEvent: "  << setw(7) << setfill('0') << fixed << nEvents;
      cout << " [" << string(barLength, '=') << ">" << "] 100.00%" << endl;
    }
  }

  float ProgressBar::getTime() {
    clock_t now = clock();
    float ret = float(now - lastTime) / CLOCKS_PER_SEC;
    lastTime = now;
    return ret;
  }

  void ProgressBar::averageTime() {
    if (nCycles < 10)
      timePerCycle = (timePerCycle * nCycles + getTime()) / (nCycles + 1);
    else
      timePerCycle = float(.98) * timePerCycle + float(.02) * getTime();
    nCycles++;
  }

  ProgressBar & ProgressBar::operator++() {

    update(++currentEvent);
    return *this;
  }
}
