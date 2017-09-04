//
// Created by mreichmann on 04.09.17.
//

#ifndef TRACKINGTELESCOPE_UTILS_H
#define TRACKINGTELESCOPE_UTILS_H

#include <string>
#include <vector>

namespace tel{

    std::string trim(const std::string &, std::string="\t\n\r\v ");
    std::vector<std::string> split(const std::string & str, const char & delim = '\t');
}

#endif //TRACKINGTELESCOPE_UTILS_H