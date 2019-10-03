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
  OutFileName(Form("ALIGNMENT/telescope%i.dat", telescope_ID)),
  PlotsDir("plots/"),
  OutDir(PlotsDir + run_number),
  FileType(".png"),
  AngleThreshold(.0005),
  ResThreshold(.000025),
  Now(clock()),
  MaximumSteps(100) {
//  MaximumSteps(14) {

  gROOT->ProcessLine("gErrorIgnoreLevel = kError;");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(53);

  fdX.resize(NPlanes, make_pair(0, 0));
  fdY.resize(NPlanes, make_pair(0, 0));
  fdA.resize(NPlanes, make_pair(0, 0));

  FR = InitFileReader();
  MaxEventNumber = FR->GetEntries();
//  MaxEventNumber = 20000;
  OrderedPlanes = GetOrderedPlanes();
  InnerPlanes = vector<unsigned short>(OrderedPlanes.begin() + 1, OrderedPlanes.end() - 1);
    std::cout<<"orderedplanes "<<std::endl;
    for(auto it=OrderedPlanes.cbegin(); it!=OrderedPlanes.cend(); it++)
        cout<<*it<<" ";
    cout<<endl;
    std::cout<<"innrplanes "<<std::endl;
    for(auto it=InnerPlanes.cbegin(); it!=InnerPlanes.cend(); it++)
        cout<<*it<<" ";
    cout<<endl;
  PlanesToAlign = AlignOnlyInnerPlanes ? InnerPlanes : vector<unsigned short>(OrderedPlanes.begin() + 1, OrderedPlanes.end());
  std::cout<<"planestoalign "<<std::endl;
  for(auto it=PlanesToAlign.cbegin(); it!=PlanesToAlign.cend(); it++)
    cout<<*it<<" ";
  cout<<endl;
  FR->GetAlignment()->ResetPlane(1, OrderedPlanes.at(0)); /** the first plane in z should be kept fixed */
  FR->GetAlignment()->SetErrorX(OrderedPlanes.at(0), 0); /** keep first plane as a fix point */
  FR->GetAlignment()->SetErrorY(OrderedPlanes.at(0), 0);
  ProgressBar = new tel::ProgressBar(MaxEventNumber);
  /** Apply Masking */
  FR->ReadPixelMask(GetMaskingFilename(TelescopeID));
  InitHistograms();

  cout << "Starting with Alignment: " << endl;
  PrintAlignment();
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
  tmp->GetAlignment()->SetErrors(TelescopeID, true);
  FILE * f = fopen("MyGainCal.dat", "w");
  tmp->GetGainCal()->PrintGainCal(f);
  fclose(f);
  return tmp;
}

void Alignment::EventLoop(const std::vector<unsigned short> & planes) {

  ProgressBar->reset();
  ResetHistograms(); /** Reset residual histograms */
  for (uint32_t i_event = 0; FR->GetNextEvent() >= 0; ++i_event) {
    if (i_event >= MaxEventNumber) break;
    ++*ProgressBar; /** print progress */

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
  for (auto i_plane: planes) {
    fdX.at(i_plane) = make_pair(hResidual.at(i_plane).GetMean(1), hResidual.at(i_plane).GetRMS(1));
    fdY.at(i_plane) = make_pair(hResidual.at(i_plane).GetMean(2), hResidual.at(i_plane).GetRMS(2));

//      double angle = atan(hResidualXdY[i_plane].GetCorrelationFactor());
    TF1 fX = TF1("fX", "pol1"), fY = TF1("fY", "pol1");
    hResidualXdY[i_plane].Fit(&fX, "q");
    hResidualYdX[i_plane].Fit(&fY, "q");
    fdA.at(i_plane) = make_pair(atan(fX.GetParameter(1)), atan(fY.GetParameter(1)));
  }

} // end EventLoop

void Alignment::PreAlign() {

  /** coarsely move the two inner planes to minimise the residuals without rotations */
  cout << "\n***************************\nPART ONE - COARSE ALIGNMENT\n***************************\n" << endl;
  FR->ResetFile();
  FR->SetPlanesUnderTest(InnerPlanes); /** ignore inner planes for tracking */

  Now = clock();
  EventLoop(InnerPlanes); /** Loop over all events and fill histograms */
  cout << Form("\nLoop duration: %2.1f\n", (clock() - Now) / CLOCKS_PER_SEC) << endl;

  for (auto i_plane: InnerPlanes){
    pair<float, float> dX(fdX.at(i_plane)), dY(fdY.at(i_plane));

    /** update the alignment */
    PLTAlignment * al = FR->GetAlignment();
    float lx(al->LX(1, i_plane)), ly(al->LY(1, i_plane));
    al->AddToLX(1, i_plane, float(dX.first));
    al->AddToLY(1, i_plane, float(dY.first));
    cout << Form("Changing alignment of plane %i in x: %+1.4f -> %+1.4f and in y: %+1.4f -> %+1.4f", i_plane, lx, al->LX(1, i_plane), ly, al->LY(1, i_plane)) << endl;

    SaveHistograms(i_plane);
  }
  PrintResiduals(InnerPlanes);
}

int Alignment::Align() {

  cout << "\n*************************\nPART TWO - FINE ALIGNMENT\n*************************\n" << endl;
  for (int i_align(0); i_align < MaximumSteps; i_align++) {
    cout << "BEGIN ITERATION " << i_align + 1 << " OUT OF " << MaximumSteps << endl;
    FR->ResetFile();
    FR->SetAllPlanes();
    EventLoop(PlanesToAlign); /** Loop over all events and fill histograms */
    int alternate = i_align / int(5);

    unsigned int elems(PlanesToAlign.size());
//    for (auto i_plane: PlanesToAlign) {
    for (unsigned int it=0; it<elems; it++) {
//      pair<float, float> dX(fdX.at(i_plane)), dY(fdY.at(i_plane)), dA(fdA.at(i_plane));
        unsigned int roc = PlanesToAlign.at(it);
        pair<float, float> dX(fdX.at(roc)), dY(fdY.at(roc)), dA(fdA.at(roc));
        if((it + alternate) % 2 == 0)
//      FR->GetAlignment()->AddToLR(1, i_plane, float(dA.first)); // DA: this was ... other_angle/3.
            FR->GetAlignment()->AddToLR(1, roc, float(dA.first - dA.second)/2.0); // DA: this was ... other_angle/3.
        else {
//      FR->GetAlignment()->AddToLX(1, i_plane, float(dX.first));
            FR->GetAlignment()->AddToLX(1, roc, float(dX.first));
//      FR->GetAlignment()->AddToLY(1, i_plane, float(dY.first));
            FR->GetAlignment()->AddToLY(1, roc, float(dY.first));
        }
//      SaveHistograms(i_plane, i_align);
      SaveHistograms(roc, i_align);
    }

    PrintResiduals(PlanesToAlign);
    cout << "END ITERATION " << i_align + 1 << " OUT OF " << MaximumSteps << endl;
//    float dX_max(GetMaxResidual(fdX)), dY_max(GetMaxResidual(fdY)), dAX_max(GetMaxAverageAngle(fdA));
    float dRes_max(GetMaximumMagRes(fdX, fdY)), dAX_max(GetMaxAverageAngle(fdA));
    cout << Form("Maximum Residuals magnitude = %f and Maximum Average Angle = %f  (ResMax->%f, AngleMax->%f)", dRes_max, dAX_max, ResThreshold, AngleThreshold) << endl;

    if (dRes_max < ResThreshold and dAX_max < AngleThreshold){
      SaveAllHistograms();
      cout << "Max angle is below " << AngleThreshold << " and max residuals are below " << ResThreshold << "=> stopping alignment." << endl;
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
    hResidual.emplace_back(        Form("Residual_ROC%i", i_plane),    Form("Residual_ROC%i",    i_plane), 1600, -4, 4, 1600, -4, 4);
//    hResidualXdY.emplace_back(TH2F(Form("ResidualXdY_ROC%i", i_plane), Form("ResidualXdY_ROC%i", i_plane), 281, -2.1075, 2.1075, 801, -1.00125, 1.00125));
    hResidualXdY.emplace_back(TProfile(Form("ResidualXdY_ROC%i", i_plane), Form("ResidualXdY_ROC%i", i_plane), 281, -2.1075, 2.1075, -1, 1));
//    hResidualYdX.emplace_back(TH2F(Form("ResidualYdX_ROC%i", i_plane), Form("ResidualYdX_ROC%i", i_plane), 401, -2.005, 2.005, 800, -1.00125, 1.00125));
    hResidualYdX.emplace_back(TProfile(Form("ResidualYdX_ROC%i", i_plane), Form("ResidualYdX_ROC%i", i_plane), 401, -2.005, 2.005, -1, 1));
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

  string sub_dir = (ind != -1) ? Form("ResROC%i/", i_plane) : "";
  string suffix = (ind != -1) ? Form("_%02i", ind) : "";
  gSystem->mkdir(OutDir + "/" + sub_dir, true);
  TCanvas Can;
  Can.cd();
  // 2D Residuals
  hResidual[i_plane].SetContour(1024);
  hResidual[i_plane].Draw("colz");
  FormatHistogram(&hResidual[i_plane], "dX [cm]", 1, "dY [cm]", 1.3, -1, 1, -1, 1);
  Can.SaveAs(OutDir + Form("/%s/%s%s%s", sub_dir.c_str(), hResidual[i_plane].GetName(), suffix.c_str(), FileType.c_str()));
  // Residual X-Projection
  gStyle->SetOptStat(1111);
  auto p_x = hResidual[i_plane].ProjectionX();
  FormatHistogram(p_x, "dX [cm]", 1, "Number of Entries", 1.3);
  p_x->Draw();
  Can.SaveAs(OutDir + Form("/%s/%s_X%s%s", sub_dir.c_str(), hResidual[i_plane].GetName(), suffix.c_str(), FileType.c_str()));
  // Residual Y-Projection
  auto p_y = hResidual[i_plane].ProjectionY();
  FormatHistogram(p_y, "dY [cm]", 1, "Number of Entries", 1.3);
  p_y->Draw();
  Can.SaveAs(OutDir + Form("/%s/%s_Y%s%s", sub_dir.c_str(), hResidual[i_plane].GetName(), suffix.c_str(), FileType.c_str()));
  // 2D Residuals X/dY
//  hResidualXdY[i_plane].SetContour(1024);
    hResidualXdY[i_plane].GetXaxis()->SetRangeUser(-1, 1);
    hResidualXdY[i_plane].GetYaxis()->SetRangeUser(-0.1, 0.1);
  hResidualXdY[i_plane].Draw();
  Can.SaveAs(OutDir + Form("/%s/%s%s%s", sub_dir.c_str(), hResidualXdY[i_plane].GetName(), suffix.c_str(), FileType.c_str()));
  // 2D Residuals Y/dX
//  hResidualYdX[i_plane].SetContour(1024);
    hResidualYdX[i_plane].GetXaxis()->SetRangeUser(-1, 1);
    hResidualYdX[i_plane].GetYaxis()->SetRangeUser(-0.1, 0.1);
    hResidualYdX[i_plane].Draw();
  Can.SaveAs(OutDir + Form("/%s/%s%s%s", sub_dir.c_str(), hResidualYdX[i_plane].GetName(), suffix.c_str(), FileType.c_str()));
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
    auto index = distance(z_pos.begin(), result);
    tmp.emplace_back(index);
    z_pos.at(index) = 999;
  }
  return tmp;
}
template <typename Q>
void Alignment::FormatHistogram(Q * h, const string & x_tit, float x_off, const string & y_tit, float y_off, float x_min, float x_max, float y_min, float y_max) {

  auto x_axis = h->GetXaxis();
  x_axis->SetTitle(x_tit.c_str());
  x_axis->SetTitleOffset(x_off);
  if (x_min != 0 and x_max != 0) { x_axis->SetRangeUser(x_min, x_max); }
  auto y_axis = h->GetYaxis();
  y_axis->SetTitle(y_tit.c_str());
  y_axis->SetTitleOffset(y_off);
  if (y_min != 0 and y_max != 0) { y_axis->SetRangeUser(y_min, y_max); }
}

void Alignment::SaveAllHistograms(int ind) {

  for (auto i_plane: OrderedPlanes){
    SaveHistograms(i_plane, ind);
  }
}

void Alignment::PrintResiduals(const vector<unsigned short> & planes) {

  cout << "\nRESIDUAL:  X         Y         aX         aY" << endl;
  for (auto i_plane: planes) {
    cout << " PLANE " << i_plane << ": " << Form("%+1.6f %+1.6f %+1.6f %+1.6f", fdX.at(i_plane).first, fdY.at(i_plane).first, fdA.at(i_plane).first, fdA.at(i_plane).second) << endl;
  }
}

float Alignment::GetMaxAverageAngle(const vector<pair<float, float>> &residuals) {
  vector<float> tmp;
  for (auto i_res: residuals){
    tmp.emplace_back(fabs((i_res.first-i_res.second))/2.0);
  }
  return *max_element(tmp.begin(), tmp.end());
}

float Alignment::GetMaximumMagRes(const vector<pair<float, float>> & residualsx, const vector<pair<float, float>> & residualsy){
    vector<float> tmp;
    unsigned long elems = residualsx.size();
    for (unsigned int i=0; i<elems; i++){
        tmp.emplace_back(sqrt(residualsx.at(i).first * residualsx.at(i).first + residualsy.at(i).first * residualsy.at(i).first));
    }
    return *max_element(tmp.begin(), tmp.end());
}
