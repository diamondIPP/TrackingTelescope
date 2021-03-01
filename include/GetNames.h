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
  static std::string type_;
  static std::vector<float> dia_z_pos_;
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
std::string GetDir();
std::string GetPlotDir();
void ValidateDirectories(const std::string&);
std::string GetCalibrationPath();
std::string GetAlignmentFilename();
std::string GetMaskingFilename();
const char * GetSignalBranchName();
uint16_t GetNPlanes();

int16_t GetRawID();
int GetNSignals();
bool UseGainInterpolator();
bool UseExternalCalibrationFunction();
bool IsROOTFile(const std::string& filename);
bool UseFileWriter();
bool FillSignalHistos();
bool UseDigitalCalibration();
template <typename T>
bool in(T , std::vector<T>);

#endif // GetNames_h
