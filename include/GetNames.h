#ifndef GetNames_h
#define GetNames_h

#include <string>

#pragma once

namespace tel {

class Config {
public:
  static int16_t telescope_id_;
  static uint16_t mask_;
  static uint16_t n_rocs_;
  static uint16_t calibration_;
  static uint16_t year_;
  static std::vector<float> dia_z_pos_;
  static int Read(int16_t);

private:
  static std::vector<float> GetZPos(uint16_t);
};

} // end tel namespace

std::string GetCalibrationFilename(int telescopeID);
std::string GetAlignmentFilename();
std::string GetDir();
uint16_t GetMaxTel();
int16_t GetRawAlignTel(uint16_t);
std::string GetMaskingFilename(int telescopeID);
const char * GetSignalBranchName();
uint8_t GetNumberOfROCS(int16_t);
int GetNumberOfSignals(int16_t);
bool GetUseGainInterpolator(int telescopeID);
bool GetUseExternalCalibrationFunction(int telescopeID);
bool GetUseRootInput(int telescopeID);
bool UseFileWriter(uint8_t);
bool FillSignalHistos(uint8_t);
bool UseDigitalCalibration(int16_t);
template <typename T>
bool in(T , std::vector<T>);

#endif // GetNames_h
