//
// Created by micha on 18.04.19.
//

#ifndef TRACKINGTELESCOPE_FINDPLANEERRORS_H
#define TRACKINGTELESCOPE_FINDPLANEERRORS_H

#include <string.h>
#include <TString.h>
class PSIFileReader;
class TF1;
class TH1F;
namespace tel { class ProgressBar; }

class FindPlaneErrors {

public:
  FindPlaneErrors(std::string, const TString &, short);
  ~FindPlaneErrors();
  short const TelescopeID;
  unsigned const NPlanes;
  float Threshold;  /** in [cm] */
  int Run();
  void SaveErrors();

private:
  std::string InFileName;
  std::string OutFileName;
  TString const PlotsDir;
  TString const OutDir;
  std::string FileType;
  PSIFileReader * FR;
  PSIFileReader * InitFileReader();
  unsigned MaxEventNumber;
  tel::ProgressBar * ProgressBar;
  unsigned short MaxIterations;
  /** Histograms */
  std::pair<TH1F*, TH1F*> hChi2All;
  std::pair<TH1F*, TH1F*> hChi2Res;
  std::vector<std::pair<float, float>> Chi2Res;
  std::pair<float, float> Chi2All;
  void FillAllChi2();
  void FillResChi2(unsigned short);
  /** Fits */
  std::pair<TF1*, TF1*> AllFit;
  std::pair<TF1*, TF1*> ResFit;
  std::pair<TF1*, TF1*> InitFits(bool=false);
  void FitGammaDistRes();
  void FitGammaDistAll();
  void SavePlots(unsigned short);
//  void SetErrors();
  void SetErrors(unsigned short);

};


#endif //TRACKINGTELESCOPE_FINDPLANEERRORS_H
