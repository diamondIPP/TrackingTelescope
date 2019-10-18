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

Alignment::Alignment(string in_file_name, const TString & run_number, short telescope_ID, bool onlyTel, unsigned short maxSteps, float maxRes, float maxAngle, unsigned long maxEvents, int8_t silDutRoc):
    TelescopeID(telescope_ID),
    NPlanes(GetNumberOfROCS(telescope_ID)),
    AlignOnlyTelescope(onlyTel),
    InFileName(std::move(in_file_name)),
    OutFileName(Form("ALIGNMENT/telescope%i.dat", telescope_ID)),
    PlotsDir("plots/"),
    OutDir(PlotsDir + run_number),
    FileType(".png"),
    AngleThreshold(maxAngle),
    ResThreshold(maxRes),
    MaxEvents(maxEvents),
    Now(clock()),
    MaximumSteps(maxSteps),
    silDUTRoc(silDutRoc){
    trackOnlyTelescope = onlyTel;

    gROOT->ProcessLine("gErrorIgnoreLevel = kError;");
    gStyle->SetOptStat(0);
    gStyle->SetPalette(53);

    fdX.resize(NPlanes, make_pair(0, 0));
    fdY.resize(NPlanes, make_pair(0, 0));
    fdA.resize(NPlanes, make_pair(0, 0));
    alignmentFinished = false;
    alignStep = 0;
    while(not alignmentFinished) {
        FR = InitFileReader();
        maxResiduals.clear();
        maxAngles.clear();
        MaxEventNumber = (MaxEvents == 0 or MaxEvents > FR->GetEntries()) ? (unsigned long) FR->GetEntries() : MaxEvents;
        OrderedPlanes = GetOrderedPlanes();
        TelescopePlanes = GetTelescopePlanes();
        DiaPlanes = GetDiamondPlanes();
        InnerPlanes = vector<unsigned short>(OrderedPlanes.begin() + 1, OrderedPlanes.end() - 1);
        std::cout << "Ordered planes: ";
        for (auto it = OrderedPlanes.cbegin(); it != OrderedPlanes.cend(); it++)
            cout << *it << " ";
        cout << endl;
        std::cout << "Telescope planes: ";
        for (auto it = TelescopePlanes.cbegin(); it != TelescopePlanes.cend(); it++)
            cout << *it << " ";
        cout << endl;
        std::cout << "Inner planes: ";
        for (auto it = InnerPlanes.cbegin(); it != InnerPlanes.cend(); it++)
            cout << *it << " ";
        cout << endl;
        std::cout << "Diamond planes: ";
        for (auto it = DiaPlanes.cbegin(); it != DiaPlanes.cend(); it++)
            cout << *it << " ";
        cout << endl;
        if(silDUTRoc != -1)
            cout << "SilDut: " << silDUTRoc << endl;

        if(alignStep == 0) {
            cout << "\n**************************************************\nAlign last telescope plane \"0\" Part 1\n**************************************************\n" << endl;
            PlanesToAlign = vector<unsigned short>(TelescopePlanes.end() - 1, TelescopePlanes.end());
            PlanesUnderTest = vector<unsigned short>(TelescopePlanes.begin() + 1, TelescopePlanes.end());
        } else if(alignStep == 1) {
            cout << "\n**************************************************\nAlign last telescope plane \"0\" Part 2 (w. tracking)\n**************************************************\n" << endl;
            PlanesToAlign = vector<unsigned short>(TelescopePlanes.end() - 1, TelescopePlanes.end());
            PlanesUnderTest = vector<unsigned short>(TelescopePlanes.begin() + 1, TelescopePlanes.end() - 1);
        } else if(alignStep == 2){
            cout << "\n**************************************************\nAlign telescope inner planes Part 2 (w. tracking)\n**************************************************\n" << endl;
            PlanesToAlign = vector<unsigned short>(TelescopePlanes.begin() + 1, TelescopePlanes.end() - 1);
            PlanesUnderTest.clear();
        } else if(alignStep == 3){
            cout << "\n**************************************************\nAlign last DUT Part 2 (w. tracking)\n**************************************************\n" << endl;
            PlanesToAlign.clear();
            PlanesToAlign.push_back((unsigned short)silDUTRoc);
            PlanesUnderTest = vector<unsigned short>(DiaPlanes.begin(), DiaPlanes.end());
        } else if(alignStep == 4) {
            cout << "\n**************************************************\nAlign DUTs\n**************************************************\n" << endl;
            PlanesToAlign = vector<unsigned short>(DiaPlanes.begin(), DiaPlanes.end());
            PlanesUnderTest = vector<unsigned short>(DiaPlanes.begin(), DiaPlanes.end());
        } else{
            cout << "Everything is aligned :)" << endl;
            break;
        }

        std::cout << "Planes to align: ";
        for (auto it = PlanesToAlign.cbegin(); it != PlanesToAlign.cend(); it++)
            cout << *it << " ";
        cout << endl;
        std::cout << "Planes under test: ";
        for (auto it = PlanesUnderTest.cbegin(); it != PlanesUnderTest.cend(); it++)
            cout << *it << " ";
        cout << endl;
        FR->GetAlignment()->ResetPlane(1, OrderedPlanes.at(0)); /** the first plane in z should be kept fixed */
        FR->GetAlignment()->SetErrorX(OrderedPlanes.at(0), 0); /** keep first plane as a fix point */
        FR->GetAlignment()->SetErrorY(OrderedPlanes.at(0), 0);
        ProgressBar = new tel::ProgressBar(MaxEventNumber);
        /** Apply Masking */
        FR->ReadPixelMask(GetMaskingFilename(TelescopeID));
        InitHistograms();

        cout << "Starting with Alignment: " << endl;
        PrintAlignment();
//        PreAlign();
        Align();
        SetNextAlignmentStep();
    }
    cout << "\nSaved plots to: " << OutDir << endl;
}

void Alignment::SetNextAlignmentStep() {
    alignStep++;
    if(alignStep == 2 and silDUTRoc == -1)
        alignStep++;
    alignmentFinished = (alignStep == 4) or (alignStep == 2 and AlignOnlyTelescope);
    FR->CloseFile();
    delete FR;
}

PSIFileReader * Alignment::InitFileReader() {
  PSIFileReader * tmp;
    if (GetUseRootInput(TelescopeID)){
        tmp = new PSIRootFileReader(InFileName, GetCalibrationFilename(TelescopeID), GetAlignmentFilename(TelescopeID), NPlanes, GetUseGainInterpolator(TelescopeID),
      GetUseExternalCalibrationFunction(TelescopeID), false, uint8_t(TelescopeID), trackOnlyTelescope);
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
      if (fabs(dR.first) >= 0.5 or fabs(dR.second) >= 0.5) continue; /** only proceed if the residual is smaller than 10mm in x or y */
      if (fabs(pow(dR.first, 2) + pow(dR.second, 2)) >= pow(0.5,2)) continue; /** only proceed if the residual is smaller than 5mm */

      if(fabs(pow((0 - dR.first) * fdY.at(i_plane).second, 2) + pow((0 - dR.second) * fdX.at(i_plane).second, 2)) > pow(nSigma * fdX.at(i_plane).second * fdY.at(i_plane).second, 2)) continue;

      hResidual[i_plane].Fill(dR.first, dR.second); // dX vs dY
      hResidualXdY[i_plane].Fill(Cluster->LX(), dR.second); // X vs dY
      hResidualYdX[i_plane].Fill(Cluster->LY(), dR.first); // Y vs dX
    }
  }
  for (auto i_plane: planes) {
    fdX.at(i_plane) = make_pair(hResidual.at(i_plane).GetMean(1), hResidual.at(i_plane).GetRMS(1));
    fdY.at(i_plane) = make_pair(hResidual.at(i_plane).GetMean(2), hResidual.at(i_plane).GetRMS(2));

      TF1 fX = TF1("fX", "pol1"), fY = TF1("fY", "pol1");
      hResidualXdY[i_plane].Fit(&fX, "q");
      hResidualYdX[i_plane].Fit(&fY, "q");
      fdA.at(i_plane) = make_pair(atan(fX.GetParameter(1)), atan(fY.GetParameter(1)));
      cout << "Plane " << i_plane << " has X: " << fdX.at(i_plane).first << " +/- " << fdX.at(i_plane).second << "; Y: " << fdY.at(i_plane).first << " +/- " << fdY.at(i_plane).second << "; Angles: X " << fdA.at(i_plane).first << "; Y " << fdA.at(i_plane).second << endl;
  }
} // end EventLoop

int Alignment::Align() {

  for (int i_align(0); i_align < MaximumSteps; i_align++) {
    cout << "BEGIN ITERATION " << i_align + 1 << " OUT OF " << MaximumSteps << endl;
    FR->ResetFile();

    FR->SetAllPlanes();
      if(not PlanesUnderTest.empty())
        FR->SetPlanesUnderTest(PlanesUnderTest);
//  nSigma shrinks with each iteration until only an ellipse of "3sigma" is used to exclude residual outliers.
      nSigma = float(3. + 1. / (1. + exp((i_align - (MaximumSteps / 7.)) / (MaximumSteps / 25.))));
      Now = clock();
    EventLoop(PlanesToAlign); /** Loop over all events and fill histograms */
      cout << Form("\nLoop duration: %2.1f\n", (clock() - Now) / CLOCKS_PER_SEC) << endl;

    unsigned int elems(PlanesToAlign.size());
    for (unsigned int it=0; it<elems; it++) {
        unsigned int roc = PlanesToAlign.at(it);
        pair<float, float> dX(fdX.at(roc)), dY(fdY.at(roc)), dA(fdA.at(roc));
            FR->GetAlignment()->AddToLR(1, roc, float(dA.first - dA.second)/2.0); // DA: calculate an average value between both corrections in Xdy and Ydx
            FR->GetAlignment()->AddToLX(1, roc, float(dX.first));
            FR->GetAlignment()->AddToLY(1, roc, float(dY.first));
      SaveHistograms(roc, i_align);
    }

    PrintResiduals(PlanesToAlign);
    cout << "END ITERATION " << i_align + 1 << " OUT OF " << MaximumSteps << endl;
    float dRes_max(GetMaximumMagRes(fdX, fdY, PlanesToAlign)), dAX_max(GetMaxAngle(fdA, PlanesToAlign));
    cout << Form("Maximum Residuals magnitude = %f and Maximum Average Angle = %f  (ResMax->%f, AngleMax->%f)", dRes_max, dAX_max, ResThreshold, AngleThreshold) << endl;
      maxResiduals.push_back(dRes_max);
      maxAngles.push_back(dAX_max);

    if (dRes_max < ResThreshold and dAX_max < AngleThreshold){
      SaveAllHistograms();
      cout << "Max residual is below " << ResThreshold << " and max angle is below " << AngleThreshold << " => stopping alignment." << endl;
      break;
    }
      if(maxResiduals.size() > 1 and maxAngles.size() > 1){
          float deltaResMax(fabs(maxResiduals.back() - maxResiduals.at(maxResiduals.size() - 2)));
          float deltaAngleMax(fabs(maxAngles.back() - maxAngles.at(maxAngles.size() - 2)));
          cout << Form("Maximum Residuals magnitude Delta = %f and Maximum Average Angle Delta = %f  (ResMaxDelta->%f, AngleMaxDelta->%f)", deltaResMax, deltaAngleMax, ResThreshold * 0.1, AngleThreshold * 0.1) << endl;
          if((deltaResMax < 0.1 * ResThreshold) and (deltaAngleMax < 0.1 * AngleThreshold)){
              SaveAllHistograms();
              cout << "Residual correction is below " << ResThreshold * 0.1 << " and max angle correction is below " << 0.1 * AngleThreshold << " => stopping alignment." << endl;
              break;
          }
      }

      cout << endl;
  } // end alignment loop

  cout << "Saving alignment file \"" << OutFileName << "\" with the following parameters:" << endl;
  PrintAlignment();
  return 0;
}

void Alignment::InitHistograms() {
  for (uint8_t i_plane(0); i_plane != NPlanes; ++i_plane){
    hResidual.emplace_back(TH2F(Form("Residual_ROC%i", i_plane),    Form("Residual_ROC%i",    i_plane), 801, -0.500625, 0.500625, 801, -0.500625, 0.500625));
    hResidualXdY.emplace_back(TProfile(Form("ResidualXdY_ROC%i", i_plane), Form("ResidualXdY_ROC%i", i_plane), 135, -0.50625, 0.50625, -1, 1));
    hResidualYdX.emplace_back(TProfile(Form("ResidualYdX_ROC%i", i_plane), Form("ResidualYdX_ROC%i", i_plane), 201, -0.5025, 0.5025, -1, 1));
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
    Can.SetGridx(); Can.SetGridy();
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
    hResidualXdY[i_plane].GetXaxis()->SetRangeUser(-0.5, 0.5);
    hResidualXdY[i_plane].GetYaxis()->SetRangeUser(-0.1, 0.1);
  hResidualXdY[i_plane].Draw();
  Can.SetGridx(); Can.SetGridy();
  Can.SaveAs(OutDir + Form("/%s/%s%s%s", sub_dir.c_str(), hResidualXdY[i_plane].GetName(), suffix.c_str(), FileType.c_str()));
  // 2D Residuals Y/dX
//  hResidualYdX[i_plane].SetContour(1024);
  hResidualYdX[i_plane].GetXaxis()->SetRangeUser(-0.5, 0.5);
  hResidualYdX[i_plane].GetYaxis()->SetRangeUser(-0.1, 0.1);
  hResidualYdX[i_plane].Draw();
  Can.SetGridx(); Can.SetGridy();
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

std::vector<unsigned short> Alignment::GetTelescopePlanes(){
    vector<unsigned short> tmp;
    vector<unsigned short> tmpfront = vector<unsigned short>(OrderedPlanes.begin(), OrderedPlanes.begin() + 2);
    vector<unsigned short> tmpback = vector<unsigned short>(OrderedPlanes.end() - 2, OrderedPlanes.end());
    tmp.reserve(tmpfront.size() + tmpback.size());
    tmp.insert(tmp.end(), tmpfront.begin(), tmpfront.end());
    tmp.insert(tmp.end(), tmpback.begin(), tmpback.end());
    return tmp;
}

std::vector<unsigned short> Alignment::GetDiamondPlanes() {
    vector<unsigned short> tmp;
    remove_copy_if(OrderedPlanes.begin(), OrderedPlanes.end(), back_inserter(tmp), [&TelescopePlanes](const unsigned short& arg){return (find(TelescopePlanes.begin(), TelescopePlanes.end(), arg) != TelescopePlanes.end());});
    if(silDUTRoc == -1)
        tmp.erase(std::remove(tmp.begin(), tmp.end(), silDUTRoc), tmp.end());
    return tmp;
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

// This method is necessary to find the plane with the maximum "angle deviation" for converging criteria
float Alignment::GetMaxAngle(const vector<pair<float, float>> &residuals, std::vector<unsigned short> alignplanes) {
  vector<float> tmp;
    unsigned long elems = residuals.size();
    tmp.resize(elems);
    for(unsigned int i = 0; i<elems; i++){
        if(std::find(alignplanes.begin(), alignplanes.end(), i) != alignplanes.end()){
//            tmp[i] = float(fabs((residuals.at(i).first - residuals.at(i).second)) / 2.0);
            tmp[i] = fabs((residuals.at(i).first));
        }
        else{
            tmp[i] = 0;
        }
    }
  return *max_element(tmp.begin(), tmp.end());
}

void Alignment::ClearVectors(){
    for(unsigned i(0); i<fdX.size(); i++){
        fdX[i] = make_pair(0, 0);
        fdY[i] = make_pair(0, 0);
        fdA[i] = make_pair(0, 0);
    }
}

// This method is necessary to find the plane with the maximum "residual magnitude" for converging criteria
float Alignment::GetMaximumMagRes(const vector<pair<float, float>> & residualsx, const vector<pair<float, float>> & residualsy, vector<unsigned short> alignplanes){
    vector<float> tmp;
    unsigned long elems = residualsx.size();
    tmp.resize(elems);
    for (unsigned int i=0; i<elems; i++){
        if(std::find(alignplanes.begin(), alignplanes.end(), i) != alignplanes.end()) {
            tmp[i] = sqrt(residualsx.at(i).first * residualsx.at(i).first + residualsy.at(i).first * residualsy.at(i).first);
        }
        else{
            tmp[i] = 0;
        }
    }
    return *max_element(tmp.begin(), tmp.end());
}
