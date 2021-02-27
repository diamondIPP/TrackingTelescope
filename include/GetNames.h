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
  static std::string plot_dir_;
  static int Read(int16_t);

private:
  static std::vector<float> GetZPos(uint16_t);
};

} // end tel namespace

struct AlignSettings {
  uint32_t max_events_ = 0;
  uint16_t n_iterations_ = 0;
  float res_thresh_ = 0;
  float angle_thresh_ = 0;
  int16_t sil_roc_ = -1;
};

AlignSettings ReadAlignSettings(std::vector<std::string>, uint16_t n_actions);
void ValidateDirectories(const std::string&);
std::string GetCalibrationFilename(int telescopeID);
std::string GetAlignmentFilename();
std::string GetDir();
std::string GetMaskingFilename();
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
