#ifndef GetNames_h
#define GetNames_h

#include <iostream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <algorithm>

std::string GetCalibrationFilename(int telescopeID);
std::string GetAlignmentFilename(int telescopeID, bool useInitial=0);
std::string GetMaskingFilename(int telescopeID);
const char * GetSignalBranchName();
int GetNumberOfROCS(int telescopeID);
int GetNumberOfSignals(int);
int GetUseGainInterpolator(int telescopeID);
int GetUseExternalCalibrationFunction(int telescopeID);
int GetUseRootInput(int telescopeID);
bool UseFileWriter(uint8_t);
bool FillSignalHistos(uint8_t);
bool UseDigitalCalibration(uint8_t);
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);
bool in(uint8_t, std::vector<uint8_t>);

#endif // GetNames_h
