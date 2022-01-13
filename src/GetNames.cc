#include "Utils.h"

#include "GetNames.h"
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <fstream>
#include "TString.h"
#include "TSystem.h"
#include "TFile.h"

using namespace std;

string GetDir() {
  /** @returns: the main directory of this program */
  string d = __FILE__;
  d = string(d, 0, d.find_last_of('/' ));
  return string(d, 0, d.find_last_of('/' ) + 1);
}

string GetPlotDir() {
  return GetDir() + "plots/";
}

namespace tel {

int16_t Config::telescope_id_;
uint16_t Config::n_rocs_, Config::mask_, Config::calibration_, Config::year_;
vector<float> Config::dia_z_pos_;
string Config::type_;

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
      s >> Config::n_rocs_ >> Config::mask_ >> Config::calibration_ >> zpos_number >> Config::year_ >> Config::type_;
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
  if (args.size() > i) {AS.max_events_ = stoi(args.at(i++));   cout << Form("Using %i events", AS.max_events_) << endl;}
  if (args.size() > i) {AS.n_iterations_ = stoi(args.at(i++)); cout << Form("Using %i iterations", AS.n_iterations_) << endl;}
  if (args.size() > i) {AS.res_thresh_ = stof(args.at(i++));   cout << Form("Using %1.2e as residual threshold", AS.res_thresh_) << endl;}
  if (args.size() > i) {AS.angle_thresh_ = stof(args.at(i++)); cout << Form("Using %1.2e as angle threshold", AS.angle_thresh_) << endl;}
  if (args.size() > i) {AS.sil_roc_ = stoi(args.at(i++));      cout << Form("Using plane %i as SIL DUT", AS.sil_roc_) << endl;}

  return AS;
}

string GetMaskingFilename(){
  /** @returns: path to the outer pixel mask file */
  string path = GetDir() + Form("data/outer_pixel_masks/%i.txt", tel::Config::mask_);
  if (gSystem->AccessPathName(path.c_str())) {
    tel::critical(Form("The mask file \"%s\" does not exist!", tel::split(path, '/').back().c_str()));
    throw;
  }
  return path;
}

string GetCalibrationPath(){
  /** @returns: path to the calibration directory */
  string path = GetDir() + Form("data/calibrations/telescope%i/", tel::Config::calibration_);
  if (gSystem->OpenDirectory(path.c_str()) == nullptr) {
    tel::critical(Form("The calibration path \"%s\" does not exist!", tel::split(path, '/').back().c_str()));
    throw;
  }
  return path + '/';
}

uint16_t GetNPlanes(){
  /** @returns: the number of the planes specified in the telescope config */
  return tel::Config::n_rocs_;
}

uint8_t GetNDUTs(){
  /** @returns: the number DUTs */
  return UseDigitalCalibration() ? size_t(GetNPlanes() - tel::Config::n_tel_planes_) : tel::Config::dia_z_pos_.size();
}

bool UseExternalCalibrationFunction() {
  /** @returns: whether to use the function given in the calibration files or not */
  return UseFileWriter();
}

bool IsROOTFile(const string & filename) {
  /** @returns: whether the file is a ROOT file or not */
  TFile f(filename.c_str());
  bool is_root_file = not f.IsZombie();
  f.Close();
  return is_root_file;
}

bool UseFileWriter(){
  /** @returns: whether to write the a root file or not */
  const int first_psi_tel = 5;
  return tel::Config::telescope_id_ >= first_psi_tel;
}

bool FillSignalHistos(){
  /** @returns: whether to fill the signal histos or not */
  const vector<int16_t> ids = {7, 8, 9};
  return in(tel::Config::telescope_id_, ids);
}

bool UseDigitalCalibration() {
  /** @returns: whether there is a pixel DUT or not */
  string type = tel::Config::type_;
  return type.find("PIX") != string::npos or type.find("Pix") != string::npos or type.find("pix") != string::npos;
}

int16_t GetRawID() {
  /** @returns: the raw telescope id, which contains only the z positions */
  const int16_t pl6(6), pl7(7), year_of_change(2016);
  if (UseDigitalCalibration()) {
    if      (tel::Config::n_rocs_ == pl6) { return -3; }
    else if (tel::Config::n_rocs_ == pl7) { return -4; }
    else    { tel::critical(Form("There is pixel raw alignment for %i planes!", tel::Config::n_rocs_)); throw; }
  } else {
    return tel::Config::year_ < year_of_change ? -1 : (tel::Config::year_ > 2020 ? -5 : -2);
  }
}

int GetNSignals() {
  /** @returns: the number of signals to show in the signal distribution. */
  return UseDigitalCalibration() ? 0 : GetNPlanes();
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

void ValidateDirectories(const string & run_number) {
  /** check if all required directories exist and create them if not */
  TString plot_dir = GetPlotDir().c_str();
  if (gSystem->OpenDirectory(plot_dir) == nullptr) {
    gSystem->mkdir(plot_dir);
  }
  if (gSystem->OpenDirectory(plot_dir + "/" + run_number) == nullptr) {
    gSystem->mkdir(plot_dir + "/" + run_number);
  }
}

bool UseGainInterpolator() { return false; }
