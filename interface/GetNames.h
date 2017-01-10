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
uint8_t GetNumberOfROCS(int16_t);
int GetNumberOfSignals(int16_t);
int GetUseGainInterpolator(int telescopeID);
int GetUseExternalCalibrationFunction(int telescopeID);
int GetUseRootInput(int telescopeID);
float GetDiamondZPosition(int16_t, uint8_t);
bool UseFileWriter(uint8_t);
bool FillSignalHistos(uint8_t);
bool UseDigitalCalibration(uint8_t);
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);
template <typename T>
bool in(T , std::vector<T>);

bool GetUseSlopeInsteadOfAngle(int16_t telescopeID);

#endif // GetNames_h
