#ifndef GetNames_h
#define GetNames_h

#include <string>

std::string GetCalibrationFilename(int telescopeID);
std::string GetAlignmentFilename();
std::string GetDir();
std::string GetMaskingFilename(int telescopeID);
const char * GetSignalBranchName();
uint8_t GetNumberOfROCS(int16_t);
int GetNumberOfSignals(int16_t);
bool GetUseGainInterpolator(int telescopeID);
bool GetUseExternalCalibrationFunction(int telescopeID);
bool GetUseRootInput(int telescopeID);
float GetDiamondZPosition(int16_t, uint8_t);
bool UseFileWriter(uint8_t);
bool FillSignalHistos(uint8_t);
bool UseDigitalCalibration(int16_t);
template <typename T>
bool in(T , std::vector<T>);

#endif // GetNames_h
