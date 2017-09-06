//
// Created by mreichmann on 04.09.17.
//

#include "Utils.h"
#include <sstream>
#include <cmath>

using namespace std;

namespace tel {

  vector<string> split(const string & s, const char & deliminator) {

    stringstream ss(s);
    string item;
    vector<string> elements;
    while (getline(ss, item, deliminator) != nullptr)
      elements.push_back(item);
    return elements;
  }

  string trim(const string &str, string trim_characters) {

    size_t b = str.find_first_not_of(trim_characters);
    size_t e = str.find_last_not_of(trim_characters);
    if (b == std::string::npos || e == std::string::npos)
      return "";
    return string(str, b, e - b + 1);
  }

  double distance(std::pair<float, float> p1, std::pair<float, float> p2){

    return sqrt(pow(p1.first - p2.first, 2) + pow(p1.second - p2.second, 2));
  }
}