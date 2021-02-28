//
// Created by micha on 18.04.19.
//

#ifndef TRACKINGTELESCOPE_FINDPLANEERRORS_H
#define TRACKINGTELESCOPE_FINDPLANEERRORS_H

#include "Action.h"
class TF1;
class TH1F;
namespace tel { class ProgressBar; }

class FindPlaneErrors : public Action {

public:
  FindPlaneErrors(const std::string & in_file_name, const TString & run_number, int16_t telescope_ID);
  ~FindPlaneErrors();
  short const TelescopeID;
  unsigned const NPlanes;
  float Threshold;  /** in [cm] */
  int Run();
  void SaveErrors();
  void PrintErrors();
  void PrintChi2s();
  void PrintChi2s(unsigned short);
  std::vector<unsigned short> GetOrderedPlanes();

private:
  std::string OutFileName;
  TString const PlotsDir;
  TString const OutDir;
  std::string FileType;
  std::vector<unsigned short> OrderedPlanes;
  unsigned MaxEventNumber;
  tel::ProgressBar * ProgressBar;
  unsigned short MaxIterations;
  /** Histograms */
  std::pair<TH1F*, TH1F*> hChi2All;
  std::pair<TH1F*, TH1F*> hChi2Res;
  std::vector<std::pair<float, float>> Chi2Res;
  std::vector<std::pair<float, float>> RealChi2Res;
  std::pair<float, float> Chi2All;
  void FillAllChi2();
  void FillResChi2();
  void FillResChi2(unsigned short, bool=true);
  /** Fits */
  std::pair<TF1*, TF1*> AllFit;
  std::pair<TF1*, TF1*> ResFit;
  std::pair<TF1*, TF1*> InitFits(bool=false);
  void FitGammaDistRes();
  void FitGammaDistAll();
  void SavePlots(unsigned short);
  void SavePlots();
  void AdjustBiggestError();
  void AdjustErrors();
  void AdjustErrors(unsigned short);

};


#endif //TRACKINGTELESCOPE_FINDPLANEERRORS_H
