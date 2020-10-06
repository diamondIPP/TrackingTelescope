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
  Threshold(.01),
  InFileName(std::move(in_file_name)),
  OutFileName(Form("ALIGNMENT/telescope%i.dat", telescope_ID)),
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

  Chi2Res.assign(NPlanes, make_pair(1, 1));
  RealChi2Res.assign(NPlanes, make_pair(1, 1));

  FR = InitFileReader();
  FR->ReadPixelMask(GetMaskingFilename(TelescopeID)); /** Apply Masking */
  OrderedPlanes = GetOrderedPlanes();
  MaxEventNumber = (FR->GetEntries() > 50000) ? 10000 : unsigned(FR->GetEntries());
  ProgressBar = new tel::ProgressBar(MaxEventNumber - 1);
  InitFits();

  Run();
}

int FindPlaneErrors::Run() {

  PrintErrors();
  FillAllChi2();
  AdjustErrors();
  PrintErrors();

  for (unsigned i(0); i < 2; i++){
    cout << "\n=========== BEGINNING ITERATION " << i + 1 << "/" << MaxIterations << " ===========\n" << endl;

    for (auto i_plane: OrderedPlanes){
      FillResChi2(i_plane);
      PrintChi2s(i_plane);
      unsigned tries(0);
      while ((fabs(Chi2Res.at(i_plane).first) > Threshold or fabs(Chi2Res.at(i_plane).second) > Threshold) and tries++ < 4){
        AdjustErrors(i_plane);
        FillResChi2(i_plane);
        PrintErrors();
        PrintChi2s(i_plane);
      }
    }
//    FillAllChi2();
//    AdjustErrors();
//    PrintErrors();
//    PrintChi2s();
//    if (fabs(Chi2All.first) < Threshold and fabs(Chi2All.second) < Threshold) { break; }
  }
  return 0;
}

void FindPlaneErrors::FillAllChi2() {

  FR->ResetFile();
  FR->SetAllPlanes();
  hChi2All.first->Reset();
  hChi2All.second->Reset();
  ProgressBar->reset();
  ProgressBar->setNEvents(MaxEventNumber);

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
  FitGammaDistAll();
  SavePlots();
  Chi2All = make_pair(AllFit.first->GetParameter(0) - 1, AllFit.second->GetParameter(0) - 1); /** fill the deviation from the expected value of 1*/
}

void FindPlaneErrors::FillResChi2() {

  cout << "Filling Chi2 with planes under test." << endl;
  ProgressBar->reset();
  ProgressBar->setNEvents(NPlanes * MaxEventNumber);
  Chi2Res.assign(NPlanes, make_pair(0, 0));
  for (unsigned short i_plane(0); i_plane < NPlanes; i_plane++){
    FillResChi2(i_plane, false);
  }
}

void FindPlaneErrors::FillResChi2(unsigned short i_plane, bool out) {

  if (out) {
    cout << Form("Filling Chi2 with plane %i under test.", i_plane) << endl;
    ProgressBar->reset();
    ProgressBar->setNEvents(MaxEventNumber);
  }
  FR->ResetFile();
  FR->SetPlaneUnderTest(i_plane);
  hChi2Res.first->Reset();
  hChi2Res.second->Reset();
  for (size_t i_event(0); FR->GetNextEvent() >= 0; ++i_event){
    if (i_event >= MaxEventNumber) break;
    ++*ProgressBar;
    if (FR->NTracks() != 1) { continue; }
    hChi2Res.first->Fill(FR->Track(0)->Chi2X());
    hChi2Res.second->Fill(FR->Track(0)->Chi2Y());
  }
  hChi2Res.first->Scale(1 / hChi2Res.first->Integral());
  hChi2Res.second->Scale(1 / hChi2Res.second->Integral());
  FitGammaDistRes();
  SavePlots(i_plane);
  Chi2Res.at(i_plane) = make_pair(ResFit.first->GetParameter(0) - 1, ResFit.second->GetParameter(0) - 1); /** fill the deviation from the expected value of 1*/
  RealChi2Res.at(i_plane) = make_pair(ResFit.first->GetChisquare(), ResFit.second->GetChisquare());
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

  ResFit.first->SetParameters(1, .2);
  ResFit.second->SetParameters(1, .2);
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

void FindPlaneErrors::SavePlots() {

  TCanvas c;
  c.cd();
  hChi2All.first->Draw();
  c.SaveAs(OutDir + "/AllPlaneChi2X" + FileType);
  hChi2All.second->Draw();
  c.SaveAs(OutDir + "/AllPlaneChi2Y" + FileType);
}

void FindPlaneErrors::AdjustBiggestError() {

  vector<float> x_diff, y_diff;
  for (auto i_chi2_res: Chi2Res){
    x_diff.emplace_back(fabs(i_chi2_res.first - Chi2All.first));
    y_diff.emplace_back(fabs(i_chi2_res.second - Chi2All.first));
  }
  size_t max_plane_x = distance(x_diff.begin(), max_element(x_diff.begin(), x_diff.end()));
  size_t max_plane_y = distance(y_diff.begin(), max_element(y_diff.begin(), y_diff.end()));
  /** increase the error for the biggest deviation. For small deviations the shape is basically not influenced */
  cout << "MAX DIFF " << x_diff.at(max_plane_x) << endl;
  for (unsigned short i_plane(0); i_plane < NPlanes; i_plane++) {
    if (i_plane != max_plane_x) { FR->GetAlignment()->SetErrorX(i_plane, FR->GetAlignment()->GetErrorX(i_plane) / (x_diff.at(max_plane_x) * .5 + 1)); }
    if (i_plane != max_plane_y) { FR->GetAlignment()->SetErrorY(i_plane, FR->GetAlignment()->GetErrorY(i_plane) / (y_diff.at(max_plane_y) * .5 + 1)); }
  }
}

void FindPlaneErrors::AdjustErrors() {

  pair<float, float> dChi2(Chi2All.first * .5, Chi2All.second * .5);
  for (unsigned short i_plane(0); i_plane < NPlanes; i_plane++) {
    FR->GetAlignment()->SetErrorX(i_plane, FR->GetAlignment()->GetErrorX(i_plane) / (dChi2.first + 1));
    FR->GetAlignment()->SetErrorY(i_plane, FR->GetAlignment()->GetErrorY(i_plane) / (dChi2.second + 1));
  }
}

void FindPlaneErrors::AdjustErrors(unsigned short plane_under_test) {

  pair<float, float> dChi2(Chi2Res.at(plane_under_test).first * .5, Chi2Res.at(plane_under_test).second * .5);
  for (auto i_plane: OrderedPlanes) {
    if (i_plane == plane_under_test) { continue; }
    FR->GetAlignment()->SetErrorX(i_plane, FR->GetAlignment()->GetErrorX(i_plane) / (dChi2.first + 1));
    FR->GetAlignment()->SetErrorY(i_plane, FR->GetAlignment()->GetErrorY(i_plane) / (dChi2.second + 1));
  }
//  FR->GetAlignment()->SetErrorX(plane_under_test, FR->GetAlignment()->GetErrorX(plane_under_test) * (dChi2.first + 1));
//  FR->GetAlignment()->SetErrorY(plane_under_test, FR->GetAlignment()->GetErrorY(plane_under_test) * (dChi2.second + 1));
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

void FindPlaneErrors::SaveErrors() {

  cout << "Saving alignment file \"" << OutFileName << "\" with the following parameters:\n" << endl;
  PrintErrors();
  FR->GetAlignment()->WriteAlignmentFile(OutFileName, NPlanes, true);
}

void FindPlaneErrors::PrintErrors() {

  PLTAlignment * al = FR->GetAlignment();
  cout << "ROC ErrorX  ErrorY" << endl;
  for (unsigned short i_plane(0); i_plane < NPlanes; i_plane++){
    cout << Form("%2i %7.4f %7.4f\n", i_plane, al->GetErrorX(i_plane), al->GetErrorY(i_plane));
  }
}

void FindPlaneErrors::PrintChi2s() {

  cout << "Chi2 ALL   (x/y): " << setprecision(4) << Chi2All.first << " / " << Chi2All.second << endl;
  for (unsigned short i_plane(0); i_plane < NPlanes; i_plane++){
    cout << "Chi2 RES " << i_plane << " (x/y): " << Chi2Res.at(i_plane).first << " / " << Chi2Res.at(i_plane).second << endl;
  }
}

void FindPlaneErrors::PrintChi2s(unsigned short i_plane) {
  cout << "Chi2 RES " << i_plane << " (x/y): " << Chi2Res.at(i_plane).first << " / " << Chi2Res.at(i_plane).second << endl;
}

std::vector<unsigned short> FindPlaneErrors::GetOrderedPlanes() {

  vector<unsigned short> tmp;
  vector<float> z_pos;
  for (unsigned i_plane(0); i_plane < NPlanes; i_plane++)
    z_pos.emplace_back(FR->GetAlignment()->LZ(1, i_plane));
  for (unsigned i_plane(0); i_plane < NPlanes; i_plane++){
    auto result = max_element(z_pos.begin(), z_pos.end());
    auto index = distance(z_pos.begin(), result);
    tmp.emplace_back(index);
    z_pos.at(index) = -999;
  }
  vector<unsigned short> res;
  res.emplace_back(tmp.front());
  res.emplace_back(tmp.back());
  for (auto val: vector<unsigned short>(tmp.begin() + 1, tmp.end() - 1)) { res.emplace_back(val); }
  return res;
}

FindPlaneErrors::~FindPlaneErrors(){

  SaveErrors();
  delete FR;
  gROOT->ProcessLine("gErrorIgnoreLevel = 0;");
}

