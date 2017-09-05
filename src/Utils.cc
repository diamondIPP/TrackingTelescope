//
// Created by mreichmann on 04.09.17.
//

#include "Utils.h"
#include <sstream>

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
}