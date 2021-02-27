#include "Utils.h"

#include "GetNames.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <fstream>
#include "TString.h"
#include "TSystem.h"

using namespace std;

uint16_t nTelescopes = 60;
vector<int16_t> pixelIDs = {10, 13, 15, 21, 25, 29, 30, 34, 35, 36, 51};
vector<int16_t> roc6IDs = {1, 2, 3, 8, 21, 25, 34};
vector<int16_t> roc7IDs = {10, 13, 15, 29, 30, 35, 51};
vector<int16_t> bcmPrimeIDs = {36, 42};


std::string GetDir() {

  string d = __FILE__;
  d = string(d, 0, d.find_last_of('/' ));
  return string(d, 0, d.find_last_of('/' ) + 1);
}

namespace tel {

int16_t Config::telescope_id_;
uint16_t Config::n_rocs_, Config::mask_, Config::calibration_, Config::year_;
vector<float> Config::dia_z_pos_;
string Config::plot_dir_ = GetDir() + "/plots";

int Config::Read(int16_t tel_id) {
  /** read the telescopes.txt config file */
  ifstream f(GetDir() + "config/telescopes.txt");
  int id, zpos_number;
  bool configured = false;
  for (string line; getline(f, line);) {
    istringstream s(line);
    s >> id;
    if (id == tel_id){
      Config::telescope_id_ = id;
      s >> Config::n_rocs_ >> Config::mask_ >> Config::calibration_ >> zpos_number >> Config::year_;
      configured = true;
    }
  }
  if (not configured){
    critical(Form("Could not find telescope %i in config file %s", tel_id, "config/telescopes.txt"));
    return 0;
  }
  Config::dia_z_pos_ = GetZPos(zpos_number);
  return 1;
}

vector<float> Config::GetZPos(uint16_t zpos_number) {
  /** read diamond z positions from config/z_pos.txt */
  ifstream f(GetDir() + "config/z_pos.txt");
  int n;
  vector<float> tmp;
  for (string line; getline(f, line);) {
    if (line.find('#') < 3) { continue; }
    line = string(line.begin(), line.begin() + line.find('#'));
    istringstream s(line);
    s >> n;
    if (n == zpos_number){
      auto v = split(trim(s.str()), ' ');
      for (const auto & word: vector<string>(v.begin() + 1, v.end())) { tmp.emplace_back(stof(word)); }
      return tmp;
    }
  }
  return tmp;
}

} // end tel namespace


AlignSettings ReadAlignSettings(vector<string> args, uint16_t n_actions) {
  /** read defualt alignment settings from config/align.txt and overwrite if args are given*/
  ifstream f(GetDir() + "config/align.txt");
  string line;
  getline(f, line); // skip first line with comments
  getline(f, line);
  istringstream s(line);

  struct AlignSettings AS;
  // read default
  s >> AS.max_events_ >> AS.n_iterations_ >> AS.res_thresh_ >> AS.angle_thresh_ >> AS.sil_roc_;
  // overwrite default if argument is provided
  uint16_t i = n_actions + 2;
  tel::warning(to_string(i));
  if (args.size() > i) {AS.max_events_ = stoi(args.at(i++));   cout << Form("Using %i events", AS.max_events_) << endl;}
  if (args.size() > i) {AS.n_iterations_ = stoi(args.at(i++)); cout << Form("Using %i iterations", AS.n_iterations_) << endl;}
  if (args.size() > i) {AS.res_thresh_ = stof(args.at(i++));   cout << Form("Using %1.2e as residual threshold", AS.res_thresh_) << endl;}
  if (args.size() > i) {AS.angle_thresh_ = stof(args.at(i++)); cout << Form("Using %1.2e as angle threshold", AS.angle_thresh_) << endl;}
  if (args.size() > i) {AS.sil_roc_ = stoi(args.at(i++));      cout << Form("Using plane %i as SIL DUT", AS.sil_roc_) << endl;}

  return AS;
}

string GetMaskingFilename(int telescopeID){

    if        (telescopeID == 1) {  return "outer_pixel_masks/outerPixelMask_Telescope1.txt";
    } else if (telescopeID == 2) {  return "outer_pixel_masks/outerPixelMask_Telescope2.txt";
    } else if (telescopeID == 5) {  return "outer_pixel_masks/outerPixelMask_Telescope5.txt";
    } else if (telescopeID == 6) {  return "outer_pixel_masks/outerPixelMask_Telescope6.txt";
    } else if (telescopeID == 7) {  return "outer_pixel_masks/outerPixelMask_Telescope7.txt";
    } else if (telescopeID == 8)  return "outer_pixel_masks/outerPixelMask_Telescope8.txt";
    else if (telescopeID == 9)  return "outer_pixel_masks/outerPixelMask_Telescope9.txt";
    else if (telescopeID == 21)  return "outer_pixel_masks/outerPixelMask_Telescope21.txt";
    else if (telescopeID == 69)  return "outer_pixel_masks/outerPixelMask_Telescope69.txt";
    else if (telescopeID == 70)  return "outer_pixel_masks/outerPixelMask_Telescope69.txt";
    else if (telescopeID >= 25)  return "outer_pixel_masks/outerPixelMask_Telescope21.txt";
    else if (telescopeID >= 10)  return "outer_pixel_masks/outerPixelMask_Telescope10.txt";
    else if (telescopeID == -1) return "outer_pixel_masks/outerPixelMask_Telescope5.txt";
    else {
        cout << "ERROR: No Masking file for telescopeID=" << telescopeID << endl;
        cout << "Exiting.." << endl;
        exit(0);
    }
}

string GetCalibrationFilename(int telescopeID){

    if      (telescopeID == 1)  return "calibration_lists/GKCalibrationList.txt";
    else if (telescopeID == 2)  return "calibration_lists/GKCalibrationList_Telescope2.txt";
    else if (telescopeID == 5)  return "calibration_lists/GKCalibrationList_Telescope5.txt";
    else if (telescopeID == 6)  return "calibration_lists/GKCalibrationList_Telescope6.txt";
    else if (telescopeID == 7)  return "calibration_lists/GKCalibrationList_Telescope7.txt";
    else if (telescopeID == 8)  return "calibration_lists/GKCalibrationList_Telescope8.txt";
    else if (telescopeID == 9)  return "calibration_lists/GKCalibrationList_Telescope9.txt";
    else if (telescopeID == 10)  return "calibration_lists/GKCalibrationList_Telescope10.txt";
    else if (telescopeID == 11)  return "calibration_lists/GKCalibrationList_Telescope9.txt";
    else if (telescopeID == 12)  return "calibration_lists/GKCalibrationList_Telescope12.txt";
    else if (telescopeID == 13)  return "calibration_lists/GKCalibrationList_Telescope13.txt";
    else if (telescopeID == 14)  return "calibration_lists/GKCalibrationList_Telescope12.txt";
    else if (telescopeID == 15)  return "calibration_lists/GKCalibrationList_Telescope13.txt";
    else if (telescopeID == 16)  return "calibration_lists/GKCalibrationList_Telescope7.txt";
    else if (telescopeID == 17)  return "calibration_lists/GKCalibrationList_Telescope12.txt";
    else if (telescopeID == 18)  return "calibration_lists/GKCalibrationList_Telescope12.txt";
    else if (telescopeID == 19)  return "calibration_lists/GKCalibrationList_Telescope12.txt";
    else if (telescopeID == 20)  return "calibration_lists/GKCalibrationList_Telescope12.txt";
    else if (telescopeID == 21)  return "calibration_lists/GKCalibrationList_Telescope21.txt";
    else if (telescopeID == 25)  return "calibration_lists/GKCalibrationList_Telescope25.txt";
    else if (telescopeID == 29)  return "calibration_lists/GKCalibrationList_Telescope29.txt";
    else if (telescopeID == 30)  return "calibration_lists/GKCalibrationList_Telescope29.txt";
    else if (telescopeID == 34)  return "calibration_lists/GKCalibrationList_Telescope34.txt";
    else if (telescopeID == 35)  return "calibration_lists/GKCalibrationList_Telescope35.txt";
    else if (telescopeID == 36)  return "calibration_lists/GKCalibrationList_Telescope35.txt";
    else if (telescopeID == 51)  return "calibration_lists/GKCalibrationList_Telescope51.txt";
    else if (telescopeID == 69)  return "calibration_lists/GKCalibrationList_Telescope69.txt";
    else if (telescopeID == 70)  return "calibration_lists/GKCalibrationList_Telescope70.txt";
    else if (telescopeID >= 10)  return "calibration_lists/GKCalibrationList_Telescope12.txt";
    else if (telescopeID == -1) return "calibration_lists/GKCalibrationList_Telescope5.txt";
    else {
        cout << "ERROR: No Calibration file for telescopeID=" << telescopeID << endl;
        cout << "Exiting.." << endl;
        exit(0);
    }
}

uint8_t GetNumberOfROCS(int16_t telescopeID){

    int16_t id = telescopeID;
    if (in(id, roc6IDs))
        return 6;
    else if (id == 4)
        return 2;
    else if (in(id, roc7IDs))
        return 7;
    else if ((id == -1) || (id >= 9))
        return 4;
    else {
        cout << "ERROR: Number of ROCs not defined for telescopeID=" << telescopeID << endl;
        cout << "Exiting.." << endl;
        exit(0);
    }
}

bool GetUseGainInterpolator(int telescopeID){

    return telescopeID == 2;
}

bool GetUseExternalCalibrationFunction(int telescopeID){

    return telescopeID == 7 || telescopeID == 8 || telescopeID == 9 || telescopeID == 10 || telescopeID >= 11;
}

bool GetUseRootInput(int telescopeID){

    return (telescopeID == -1) || (telescopeID == 7) || telescopeID == 8 || telescopeID == 9 || telescopeID == 10 || telescopeID >= 11;
}

bool UseFileWriter(uint8_t telescopeID){

    vector<uint8_t> ids;
    for (uint8_t id = 7; id <= nTelescopes; id++)
        ids.push_back(id);
    return in(telescopeID, ids);
}

bool FillSignalHistos(uint8_t telescopeID){

    vector<uint8_t> ids = {7, 8, 9};
    return in(telescopeID, ids);
}

bool UseDigitalCalibration(int16_t telescopeID){

    return in(telescopeID, pixelIDs);

}

int GetNumberOfSignals(int16_t telescopeID){

    int16_t id = telescopeID;
    if (in(id, pixelIDs))
        return 0;
    else if ((id == 7) || (id == 8) || (id == 9) || id >= 11)
        return 4;
    else {
        cerr << "ERROR: Number of Signals is not defined for telescopeID=" << telescopeID << endl;
        cout << "Exiting.." << endl;
        exit(0);
    }
}

const char * GetSignalBranchName(){

    return "blub";
}

template<typename T>
bool in(T num, vector<T> ids){

    return find(ids.begin(), ids.end(), num ) != ids.end();
}


/** Get the correct alignment for a given telescope */
string GetAlignmentFilename(){

  return GetDir() + "data/alignments.txt";
}

uint16_t GetMaxTel() {
  /** @returns highest used telescope number */
  vector<int16_t> numbers;
  ifstream f(GetAlignmentFilename());
  int tel;
  for (string line; getline(f, line);){
    istringstream s(line);
    s >> tel;
    numbers.emplace_back(tel);
  }
  return *max_element(numbers.begin(), numbers.end());
}

int16_t GetRawAlignTel(uint16_t n_planes) {

  if (n_planes == 4) { return -2; }
  else if (n_planes == 6) { return -3; }
  else if (n_planes == 7) { return -4; }
  else { return 0; }
}

void ValidateDirectories(const string & run_number) {
  /** check if all required directories exist and create them if not */
  TString plot_dir = tel::Config::plot_dir_.c_str();
  if (gSystem->OpenDirectory(plot_dir) == nullptr) {
    gSystem->mkdir(plot_dir);
  }
  if (gSystem->OpenDirectory(plot_dir + "/" + run_number) == nullptr) {
    gSystem->mkdir(plot_dir + "/" + run_number);
  }
}
