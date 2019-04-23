//
// Created by micha on 18.04.19.
//

#include "FindPlaneErrors.h"
#include "GetNames.h"
#include "PSIBinaryFileReader.h"
#include "PSIRootFileReader.h"
#include "TF1.h"
#include "Utils.h"
#include "TH1F.h"

using namespace std;

FindPlaneErrors::FindPlaneErrors(string in_file_name, const TString & run_number, short telescope_ID):
  TelescopeID(telescope_ID),
  NPlanes(GetNumberOfROCS(telescope_ID)),
  Threshold(.03),
  InFileName(std::move(in_file_name)),
  OutFileName("NewAlignment.dat"),
  PlotsDir("plots/"),
  OutDir(PlotsDir + run_number),
  FileType(".png"),
  MaxIterations(8),
  hChi2All(new TH1F("Chi2XAll", "Chi2XAll", 100, 0, 20), new TH1F("Chi2YAll", "Chi2YAll", 100, 0, 20)),
  hChi2Res(new TH1F("Chi2XRes", "Chi2XRes", 100, 0, 20), new TH1F("Chi2YRes", "Chi2YRes", 100, 0, 20)),
  Chi2All(-999, -999),
  AllFit(InitFits()),
  ResFit(InitFits(true)) {

  gROOT->ProcessLine("gErrorIgnoreLevel = kError;");

  FR = InitFileReader();
  FR->ReadPixelMask(GetMaskingFilename(TelescopeID)); /** Apply Masking */
  MaxEventNumber = unsigned(FR->GetEntries());
//  MaxEventNumber = 10000;
  ProgressBar = new tel::ProgressBar(MaxEventNumber - 1);
  InitFits();

  Run();
}

int FindPlaneErrors::Run() {

  for (unsigned i(0); i < MaxIterations; i++){
    cout << "\n=========== BEGINNING ITERATION " << i + 1 << "/" << MaxIterations << " ===========\n" << endl;

    FillAllChi2(); /** Fill All Plane Chi2 */
    FitGammaDistAll();
    TCanvas c;
    c.cd();
    hChi2All.first->Draw();
    c.SaveAs(OutDir + "/blub.png");
    Chi2All = make_pair(AllFit.first->GetParameter(0), AllFit.second->GetParameter(0));
    cout << "Chi2 of all planes (x/y): " << setprecision(4) << Chi2All.first << " / " << Chi2All.second << endl;

    Chi2Res.resize(NPlanes, make_pair(0, 0));
    for (unsigned short i_plane(0); i_plane < NPlanes; i_plane++){
      while ((fabs(1 - Chi2Res.at(i_plane).first) > Threshold) or (fabs(1 - Chi2Res.at(i_plane).second) > Threshold)) {
        FillResChi2(i_plane); /** Fill Chi2 with one plane under test */
        FitGammaDistRes();
        SavePlots(i_plane);
        Chi2Res.at(i_plane) = make_pair(ResFit.first->GetParameter(0), ResFit.second->GetParameter(0));
        SetErrors(i_plane);
      }
    }
    bool done(true);
    cout << "Chi2 of all planes (x/y): " << Chi2All.first << " / " << Chi2All.second << endl;
    for (unsigned short i_plane(0); i_plane < NPlanes; i_plane++){
      done = (fabs(1 - Chi2Res.at(i_plane).first) < Threshold) and (fabs(1 - Chi2Res.at(i_plane).second) < Threshold);
      cout << "Chi2 for plane " << i_plane << " (x/y): " << 1 - Chi2Res.at(i_plane).first << " / " << 1 - Chi2Res.at(i_plane).second << endl;
    }
    if (done) { break; }
  }
  return 0;
}

void FindPlaneErrors::FillAllChi2() {

  FR->ResetFile();
  FR->SetAllPlanes();
  hChi2All.first->Reset();
  hChi2All.second->Reset();

  cout << "\nFilling chi2 histograms with all planes tracked." << endl;\
  for (size_t i_event(0); FR->GetNextEvent() >= 0; ++i_event){
    if (i_event >= MaxEventNumber) break;
    if (FR->NTracks() != 1) { continue; }
    ProgressBar->update(i_event);
    hChi2All.first->Fill(FR->Track(0)->Chi2X());
    hChi2All.second->Fill(FR->Track(0)->Chi2Y());
  }
  hChi2All.first->Scale(1 / hChi2All.first->Integral());
  hChi2All.second->Scale(1 / hChi2All.second->Integral());
}

void FindPlaneErrors::FillResChi2(unsigned short i_plane) {

  ProgressBar->reset();
  FR->ResetFile();
  FR->SetPlaneUnderTest(i_plane);
  hChi2Res.first->Reset();
  hChi2Res.second->Reset();

  cout << "Filling Chi2 histogram with plane " << i_plane << " under test." << endl;
  for (size_t i_event(0); FR->GetNextEvent() >= 0; ++i_event){
    if (i_event >= MaxEventNumber) break;
    ProgressBar->update(i_event);
    if (FR->NTracks() != 1) { continue; }
    hChi2Res.first->Fill(FR->Track(0)->Chi2X());
    hChi2Res.second->Fill(FR->Track(0)->Chi2Y());
  }
  hChi2Res.first->Scale(1 / hChi2Res.first->Integral());
  hChi2Res.second->Scale(1 / hChi2Res.second->Integral());
}

pair<TF1*, TF1*> FindPlaneErrors::InitFits(bool dut) {

  unsigned short ndf = NPlanes - (dut ? 3 : 2);
  float max_chi2 = dut ? 10 : 20;
  vector<TF1*> tmp;
  for (unsigned i_hist(0); i_hist < 2; i_hist++){
    tmp.emplace_back(new TF1(Form("f%i", i_hist + (dut ? 2 : 0)), Form("[1] * TMath::GammaDist([0] * x, %i / 2., 0, 2)", ndf), 0, max_chi2));
    tmp.at(i_hist)->SetNpx(1000);
    tmp.at(i_hist)->SetParameters(1, .2);
    tmp.at(i_hist)->SetParLimits(0, .1, 10);
    tmp.at(i_hist)->SetParLimits(1, .01, 100);
  }
  return make_pair(tmp.at(0), tmp.at(1));
}

void FindPlaneErrors::FitGammaDistAll() {

  AllFit.first->SetParameters(1, .2);
  AllFit.second->SetParameters(1, .2);
  hChi2All.first->Fit(AllFit.first, "q");
  hChi2All.second->Fit(AllFit.second, "q");
}

void FindPlaneErrors::FitGammaDistRes() {

  AllFit.first->SetParameters(0, 1, 1, 1);
  AllFit.second->SetParameters(0, 1, 1, 1);
  hChi2Res.first->Fit(ResFit.first, "q");
  hChi2Res.second->Fit(ResFit.second, "q");
}

void FindPlaneErrors::SavePlots(unsigned short i_plane) {

  TCanvas c;
  c.cd();
  hChi2Res.first->SetLineColor(3);
  hChi2Res.second->SetLineColor(3);
  hChi2Res.first->Draw();
  hChi2All.first->Draw("same");
  c.SaveAs(OutDir + Form("/FunWithChi2X_ROC%i", i_plane) + FileType);
  hChi2Res.second->Draw();
  hChi2All.second->Draw("same");
  c.SaveAs(OutDir + Form("/FunWithChi2Y_ROC%i", i_plane) + FileType);
}

//void FindPlaneErrors::SetErrors() {
//
//  pair<float, float> max_dchi2(-1, -1);
//  pair<unsigned short, unsigned short> max_ind(-1, -1);
//  for (unsigned short i_plane(0); i_plane < NPlanes; i_plane++){
//    cout << "Chi2 for plane " << i_plane << ": " << Chi2Res.at(i_plane).first << " in x and " << Chi2Res.at(i_plane).second << " in y" << endl;
//    if (fabs(Chi2All.first - Chi2Res.at(i_plane).first) > max_dchi2.first){
//      max_ind.first = i_plane;
//      max_dchi2.first = fabs(Chi2All.first - Chi2Res.at(i_plane).first);
//    }
//    if (fabs(Chi2All.second - Chi2Res.at(i_plane).second) > max_dchi2.second){
//      max_ind.second = i_plane;
//      max_dchi2.second = fabs(Chi2All.second - Chi2Res.at(i_plane).second);
//    }
//  }
//  FR->GetAlignment()->SetErrorX(max_ind.first, FR->GetAlignment()->GetErrorX(max_ind.first) / fabs(Chi2All.first - Chi2Res.at(max_ind.first).first) / 2);
//  FR->GetAlignment()->SetErrorY(max_ind.second, FR->GetAlignment()->GetErrorY(max_ind.second) / fabs(Chi2All.second - Chi2Res.at(max_ind.second).second) / 2);
//
//  for (unsigned short i_plane(0); i_plane < NPlanes; i_plane++){
//    FR->GetAlignment()->SetErrorX(i_plane, FR->GetAlignment()->GetErrorX(i_plane) / Chi2All.first);
//    FR->GetAlignment()->SetErrorY(i_plane, FR->GetAlignment()->GetErrorY(i_plane) / Chi2All.second);
//    cout <<  "Residual for plane " << i_plane << ": " << FR->GetAlignment()->GetErrorX(i_plane) << " in x and " << FR->GetAlignment()->GetErrorY(i_plane) << " in y" << endl;
//  }
//}

void FindPlaneErrors::SetErrors(unsigned short i_pl) {

  pair<float, float> dChi2((Chi2Res.at(i_pl).first - 1) / (NPlanes - 1), (Chi2Res.at(i_pl).second - 1) / (NPlanes - 1));
  cout << "dChi2 for plane " << i_pl << " (x/y): " << Chi2Res.at(i_pl).first - 1 << " / " << Chi2Res.at(i_pl).second - 1 << endl;
  for (unsigned short i_plane(0); i_plane < NPlanes; i_plane++) {
    if (i_plane == i_pl) { continue; }
    FR->GetAlignment()->SetErrorX(i_plane, FR->GetAlignment()->GetErrorX(i_plane) / (dChi2.first + 1));
    FR->GetAlignment()->SetErrorY(i_plane, FR->GetAlignment()->GetErrorY(i_plane) / (dChi2.second + 1));
    cout << "Residual for plane " << i_plane << " (x/y): " << FR->GetAlignment()->GetErrorX(i_plane) << " / " << FR->GetAlignment()->GetErrorY(i_plane) << endl;
  }
}


PSIFileReader * FindPlaneErrors::InitFileReader() {

  PSIFileReader * tmp;
  if (GetUseRootInput(TelescopeID)){
    tmp = new PSIRootFileReader(InFileName, GetCalibrationFilename(TelescopeID), GetAlignmentFilename(TelescopeID), NPlanes, GetUseGainInterpolator(TelescopeID),
                                GetUseExternalCalibrationFunction(TelescopeID), false, uint8_t(TelescopeID));
  }
  else {
    tmp = new PSIBinaryFileReader(InFileName, GetCalibrationFilename(TelescopeID), GetAlignmentFilename(TelescopeID), NPlanes, GetUseGainInterpolator(TelescopeID),
                                  GetUseExternalCalibrationFunction(TelescopeID));
  }
  tmp->GetAlignment()->SetErrors(TelescopeID, true);
//  tmp->GetAlignment()->SetErrors(TelescopeID);
  FILE * f = fopen("MyGainCal.dat", "w");
  tmp->GetGainCal()->PrintGainCal(f);
  fclose(f);
  return tmp;
}

FindPlaneErrors::~FindPlaneErrors(){

  SaveErrors();
  delete FR;
  gROOT->ProcessLine("gErrorIgnoreLevel = 0;");
}

void FindPlaneErrors::SaveErrors() {

  PLTAlignment * al = FR->GetAlignment();
  cout << "Saving alignment file \"" << OutFileName << "\" with the following parameters:\n" << endl;
  cout << "ROC ErrorX  ErrorY" << endl;
  for (unsigned short i_plane(0); i_plane < NPlanes; i_plane++){
    cout << Form("%2i %7.4f %7.4f\n", i_plane, al->GetErrorX(i_plane), al->GetErrorY(i_plane));
  }
  al->WriteAlignmentFile(OutFileName, NPlanes, true);
}

