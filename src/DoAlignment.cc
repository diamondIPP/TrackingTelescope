#include <utility>

#include <Utils.h>
#include <DoAlignment.h>
#include "PSIBinaryFileReader.h"
#include "PSIRootFileReader.h"
#include "GetNames.h"
#include "TestPlaneEfficiencySilicon.h"
#include "PLTPlane.h"
#include "PLTAlignment.h"

using namespace std;

Alignment::Alignment(string in_file_name, const TString & run_number, short telescope_ID, bool align_only_inner_planes):
  TelescopeID(telescope_ID),
  NPlanes(GetNumberOfROCS(telescope_ID)),
  AlignOnlyInnerPlanes(align_only_inner_planes),
  InFileName(std::move(in_file_name)),
  OutFileName("NewAlignment.dat"),
  PlotsDir("plots/"),
  OutDir(PlotsDir + run_number),
  FileType(".png"),
  AngleThreshold(.01),
  TotResThreshold(.01),
  Now(clock()),
  MaximumSteps(14) {

  gROOT->ProcessLine("gErrorIgnoreLevel = kError;");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(53);

  FR = InitFileReader();
  MaxEventNumber = unsigned(dynamic_cast<PSIRootFileReader*>(FR)->fTree->GetEntries());
//  MaxEventNumber = 20000;
  OrderedPlanes = GetOrderedPlanes();
  InnerPlanes = vector<unsigned short>(OrderedPlanes.begin() + 1, OrderedPlanes.end() - 1);
  PlanesToAlign = AlignOnlyInnerPlanes ? InnerPlanes : vector<unsigned short>(OrderedPlanes.begin() + 1, OrderedPlanes.end());
  FR->GetAlignment()->ResetPlane(1, OrderedPlanes.at(0)); /** the first plane in z should be kept fixed */
  ProgressBar = new tel::ProgressBar(MaxEventNumber - 1);
  /** Apply Masking */
  FR->ReadPixelMask(GetMaskingFilename(TelescopeID));
  InitHistograms();
    
  PreAlign();
  Align();
  cout << "\nSaved plots to: " << OutDir << endl;
}

PSIFileReader * Alignment::InitFileReader() {

  PSIFileReader * tmp;
  if (GetUseRootInput(TelescopeID)){
    tmp = new PSIRootFileReader(InFileName, GetCalibrationFilename(TelescopeID), GetAlignmentFilename(TelescopeID), NPlanes, GetUseGainInterpolator(TelescopeID),
      GetUseExternalCalibrationFunction(TelescopeID), false, uint8_t(TelescopeID));
  }
  else {
    tmp = new PSIBinaryFileReader(InFileName, GetCalibrationFilename(TelescopeID), GetAlignmentFilename(TelescopeID), NPlanes, GetUseGainInterpolator(TelescopeID),
      GetUseExternalCalibrationFunction(TelescopeID));
  }
  tmp->GetAlignment()->SetErrors(TelescopeID);
  FILE * f = fopen("MyGainCal.dat", "w");
  tmp->GetGainCal()->PrintGainCal(f);
  fclose(f);
  return tmp;
}

void Alignment::EventLoop(const std::vector<unsigned short> & planes) {

  for (uint32_t i_event = 0; FR->GetNextEvent() >= 0; ++i_event) {
    if (i_event >= MaxEventNumber) break;
    ProgressBar->update(i_event); /** print progress */

    if (FR->NTracks() != 1) continue; /** proceed only if we have exactly one track */
    PLTTrack * Track = FR->Track(0);

    for (auto i_plane: planes) {

      PLTPlane * Plane = FR->Plane(i_plane);
      if (Plane->NClusters() != 1) continue; /** proceed only if there is exactly one cluster */

      PLTCluster * Cluster = Plane->Cluster(0);
      /** if the plane is under test, calculate the residuals otherwise just take them from the track */
      pair<float, float> dR = (FR->IsPlaneUnderTest(i_plane)) ? Track->GetResiduals(*Cluster, *FR->GetAlignment()) : Track->GetResiduals(i_plane);
      if (fabs(dR.first) >= 2 or fabs(dR.second) >= 2) continue; /** only proceed if the residual is smaller than 2mm in x or y */

      hResidual[i_plane].Fill(dR.first, dR.second); // dX vs dY
      hResidualXdY[i_plane].Fill(Cluster->LX(), dR.second); // X vs dY
      hResidualYdX[i_plane].Fill(Cluster->LY(), dR.first); // Y vs dX
    }
  }
} // end EventLoop

void Alignment::PreAlign() {

  /** coarsely move the two inner planes to minimise the residuals without rotations */
  cout << "\n***************************\nPART ONE - COARSE ALIGNMENT\n***************************\n" << endl;
  FR->ResetFile();
  FR->SetPlanesUnderTest(InnerPlanes); /** ignore inner planes for tracking */
  ResetHistograms(); /** Reset residual histograms */
  PrintAlignment();

  Now = clock();
  EventLoop(InnerPlanes); /** Loop over all events and fill histograms */
  cout << "\nLoop duration:" << (clock() - Now) / CLOCKS_PER_SEC << endl;
  cout << "RES " << hResidual.at(1).GetMean() << endl;

  for (auto i_plane: InnerPlanes){
    pair<float, float> dR = GetMeanResiduals(i_plane);
    pair<float, float> rms = GetRMS(i_plane);
    cout << "\nPLANE " << i_plane << "\n--RESIDUALS: " << setprecision(3) << dR.first << " (" << rms.first << "), " << dR.second << " (" << rms.second << ")" << endl;

    /** update the alignment */
    cout << "Before: " << FR->GetAlignment()->LX(1, i_plane) << " / " << FR->GetAlignment()->LY(1, i_plane) << endl;
    FR->GetAlignment()->AddToLX(1, i_plane, float(dR.first));
    FR->GetAlignment()->AddToLY(1, i_plane, float(dR.second));
    cout << "After:  " << FR->GetAlignment()->LX(1, i_plane) << " / " << FR->GetAlignment()->LY(1, i_plane) << endl;

    SaveHistograms(i_plane);
  }
}

int Alignment::Align() {

  cout << "\n*************************\nPART TWO - FINE ALIGNMENT\n*************************\n" << endl;
  for (int i_align(0); i_align < MaximumSteps; i_align++) {
    cout << "BEGIN ITERATION " << i_align + 1 << " OUT OF " << MaximumSteps << endl;
    FR->ResetFile();
    FR->SetAllPlanes();
    ProgressBar->reset();

    EventLoop(PlanesToAlign); /** Loop over all events and fill histograms */

    float total_angle = 0;
    float total_res = 0;

    for (auto i_plane: PlanesToAlign) {

      double x_res = hResidual[i_plane].GetMean(1), x_rms = hResidual[i_plane].GetRMS(1);
      double y_res = hResidual[i_plane].GetMean(2), y_rms = hResidual[i_plane].GetRMS(2);
      cout << "\nPLANE " << i_plane << "\n--RESIDUALS: " << setprecision(3) << x_res << " (" << x_rms << "), " << y_res << " (" << y_rms << ")" << endl;

      double angle = atan(hResidualXdY[i_plane].GetCorrelationFactor());
      TF1 fX = TF1("fX", "pol1"), fY = TF1("fY", "pol1");
      hResidualXdY[i_plane].Fit(&fX, "q");
      hResidualYdX[i_plane].Fit(&fY, "q");

      double angle_x = atan(fX.GetParameter(1)), angle_y = atan(fY.GetParameter(1));
      total_angle += fabs(angle_x);
      total_res += fabs(x_res) + fabs(y_res);
      cout << "--Angle Corr: " << angle << ", Angle XdY: " << angle_x << ", Angle YdX: " << angle_y << endl;

      SaveHistograms(i_plane, (total_angle < AngleThreshold && total_res < TotResThreshold) ? -1 : i_align);

      FR->GetAlignment()->AddToLR(1, i_plane, float(angle_x)); // DA: this was ... other_angle/3.
      FR->GetAlignment()->AddToLX(1, i_plane, float(x_res));
      FR->GetAlignment()->AddToLY(1, i_plane, float(y_res));
    }
    cout << "END ITERATION " << i_align + 1 << " OUT OF 14" << endl;
    cout << "\nSum of angle corrections:    " << total_angle << endl;
    cout << "Sum of residuals in x and y: " << total_res << endl;

    if (total_angle < AngleThreshold && total_res < TotResThreshold){
      std::cout << "Total angle is below " << AngleThreshold << " and Total Residual is below " << TotResThreshold << "=> stopping alignment." << endl;
      break;
    }
    cout << endl;
  } // end alignment loop

  cout << "Saving alignment file \"" << OutFileName << "\" with the following parameters:" << endl;
  PrintAlignment();
  return 0;
}

void Alignment::InitHistograms() {

  for (uint8_t i_plane(0); i_plane != NPlanes; ++i_plane){
    hResidual.emplace_back(        Form("Residual_ROC%i", i_plane),    Form("Residual_ROC%i",    i_plane), 400, -.5, .5, 400, -.5, .5);
    hResidualXdY.emplace_back(TH2F(Form("ResidualXdY_ROC%i", i_plane), Form("ResidualXdY_ROC%i", i_plane), 133, -1, 0.995, 400, -.5, .5));
    hResidualYdX.emplace_back(TH2F(Form("ResidualYdX_ROC%i", i_plane), Form("ResidualYdX_ROC%i", i_plane), 201, -1, 1, 400, -.5, .5));
  }
}

void Alignment::ResetHistograms() {

  for (uint8_t i_plane(0); i_plane != NPlanes; ++i_plane){
    hResidual.at(i_plane).Reset();
    hResidualXdY.at(i_plane).Reset();
    hResidualYdX.at(i_plane).Reset();
  }
}

void Alignment::SaveHistograms(unsigned i_plane, int ind) {

  TString sub_dir = (ind != -1) ? to_string(ind) + "/" : "";
  gSystem->mkdir(OutDir + "/" + sub_dir, true);
  TCanvas Can;
  Can.cd();
  // 2D Residuals
  hResidual[i_plane].SetContour(1024);
  hResidual[i_plane].Draw("colz");
  Can.SaveAs(OutDir + "/" + sub_dir + TString(hResidual[i_plane].GetName()) + FileType);
  // Residual X-Projection
  gStyle->SetOptStat(1111);
  auto p_x = hResidual[i_plane].ProjectionX();
  FormatHistogram(p_x, "dX [cm]", 1, "Number of Entries", 1.3);
  p_x->Draw();
  Can.SaveAs(OutDir + "/" + sub_dir + TString(hResidual[i_plane].GetName()) + "_X" + FileType);
  // Residual Y-Projection
  auto p_y = hResidual[i_plane].ProjectionY();
  FormatHistogram(p_y, "dY [cm]", 1, "Number of Entries", 1.3);
  p_y->Draw();
  Can.SaveAs(OutDir + "/" + sub_dir + TString(hResidual[i_plane].GetName()) + "_Y" + FileType);
  // 2D Residuals X/dY
  hResidualXdY[i_plane].SetContour(1024);
  hResidualXdY[i_plane].Draw("colz");
  Can.SaveAs(OutDir + "/" + sub_dir + TString(hResidualXdY[i_plane].GetName()) + FileType);
  // 2D Residuals Y/dX
  hResidualYdX[i_plane].SetContour(1024);
  hResidualYdX[i_plane].Draw("colz");
  Can.SaveAs(OutDir + "/" + sub_dir + TString(hResidualYdX[i_plane].GetName()) + FileType);
}

Alignment::~Alignment() {

  delete FR;
  gROOT->ProcessLine("gErrorIgnoreLevel = 0;");
}

void Alignment::PrintAlignment() {

  PLTAlignment * al = FR->GetAlignment();
  for (unsigned i_plane(0); i_plane < NPlanes; i_plane++)
    cout << Form("%2i %1i %15E %15E %15E %15E\n", 1, i_plane, al->LR(1, i_plane), al->LX(1, i_plane), al->LY(1, i_plane), al->LZ(1, i_plane));
  al->WriteAlignmentFile(OutFileName, FR->NMAXROCS);

}

std::vector<unsigned short> Alignment::GetOrderedPlanes() {

  vector<unsigned short> tmp;
  vector<float> z_pos;
  for (unsigned i_plane(0); i_plane < NPlanes; i_plane++)
    z_pos.emplace_back(FR->GetAlignment()->LZ(1, i_plane));
  for (unsigned i_plane(0); i_plane < NPlanes; i_plane++){
    auto result = min_element(z_pos.begin(), z_pos.end());
    tmp.emplace_back(distance(z_pos.begin(), result));
    z_pos.erase(result);
  }
  return tmp;
}

void Alignment::FormatHistogram(TH1 * h, const string & x_tit, float x_off, const string & y_tit, float y_off) {

  auto x_axis = h->GetXaxis();
  x_axis->SetTitle(x_tit.c_str());
  x_axis->SetTitleOffset(x_off);
  auto y_axis = h->GetYaxis();
  y_axis->SetTitle(y_tit.c_str());
  y_axis->SetTitleOffset(y_off);
}

std::pair<float, float> Alignment::GetMeanResiduals(unsigned short i_plane) { return std::make_pair(hResidual[i_plane].GetMean(1), hResidual[i_plane].GetMean(2)); }
std::pair<float, float> Alignment::GetRMS(unsigned short i_plane) { return std::make_pair(hResidual[i_plane].GetRMS(1), hResidual[i_plane].GetRMS(2)); }
