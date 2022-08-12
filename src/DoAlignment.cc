#include <utility>
#include <algorithm>

#include <Utils.h>
#include <DoAlignment.h>
#include "GetNames.h"
#include "TestPlaneEfficiencySilicon.h"
#include "PLTPlane.h"

using namespace std;

Alignment::Alignment(const string & in_file_name, const TString & run_number, uint16_t telescope_id, bool only_tel, uint16_t max_steps,
                     float max_res, float max_angle, uint32_t max_events, int16_t sil_dut_roc):
  Action(in_file_name, run_number),
  telescope_id_(telescope_id),
  n_planes_(GetNPlanes()),
  align_only_telescope_(only_tel),
  at_step_(0),
  alignment_finished_(false),
  plots_dir_(GetPlotDir() + run_number),
  file_type_(".png"),
  angle_thresh_(max_angle),
  res_thresh_(max_res),
  now_(clock()),
  maximum_steps_(max_steps),
  sil_dut_roc_(sil_dut_roc),
  max_sigma_(4), min_sigma_(3), n_sigma_(max_sigma_) {

    ConfigROOT();

    fdX.resize(n_planes_, make_pair(0, 1));
    fdY.resize(n_planes_, make_pair(0, 1));
    fdA.resize(n_planes_, make_pair(0, 1));

    /** set the correct order of planes */
    FR = InitFileReader();
    max_event_number_ = max_events == 0 ? FR->GetEntries() : min(max_events, FR->GetEntries());
    cout << endl << "Found plane configuration: " << endl;
    ordered_planes_ = GetOrderedPlanes();
    telescope_planes_ = GetTelescopePlanes();
    dia_planes_ = GetDiamondPlanes();
    inner_planes_ = vector<uint16_t>(ordered_planes_.begin() + 1, ordered_planes_.end() - 1);
    cout << "Inner planes: " << tel::to_string(inner_planes_) << endl;

    while (not alignment_finished_) {

      SetPlanes();
      if (planes_to_align_.empty()){
        cout << endl << "Nothing to do ..." << endl;
        SetNextAlignmentStep();
        continue;
      }

      FR = InitFileReader();
      last_max_res_ = make_pair(0, 0);
      delta_max_res_ = make_pair(0, 0);


      /** keep first plane as a fix point */
      FR->GetAlignment()->ResetPlane(1, ordered_planes_.at(0));
      FR->GetAlignment()->SetErrorX(ordered_planes_.at(0), 0);
      FR->GetAlignment()->SetErrorY(ordered_planes_.at(0), 0);

      ProgressBar = new tel::ProgressBar(max_event_number_);

      /** Apply Masking */
      FR->ReadPixelMask(GetMaskingFilename());
      InitHistograms();

      cout << "\nStarting with Alignment: " << endl;
      PrintAlignment();
      Align();
      SetNextAlignmentStep();

      FR->CloseFile();
      delete FR;
    }

} // end constructor

Alignment::~Alignment() {

  gROOT->ProcessLine("gErrorIgnoreLevel = 0;");
}

void Alignment::SetPlanes() {
  planes_to_align_.clear();
  planes_under_test_.clear();
  if(at_step_ == 0) {
//    tel::print_banner("PART I: Align last telescope plane");
    tel::print_banner("PART I: Translate all planes relative to first plane in beam");
//    planes_to_align_ = vector<uint16_t>(telescope_planes_.end() - 1, telescope_planes_.end());
    planes_to_align_ = vector<uint16_t>(ordered_planes_.begin() + 1, ordered_planes_.end());
    planes_under_test_ = vector<uint16_t>(telescope_planes_.begin() + 1, telescope_planes_.end());
  } else if(at_step_ == 1) {
    tel::print_banner("PART II: Align last telescope plane (w. tracking)");
    planes_to_align_ = vector<uint16_t>(telescope_planes_.end() - 1, telescope_planes_.end());
    planes_under_test_ = vector<uint16_t>(telescope_planes_.begin() + 1, telescope_planes_.end() - 1);
  } else if(at_step_ == 2){
    tel::print_banner("PART III: Align inner telescope planes (w. tracking)");
    planes_to_align_ = vector<uint16_t>(telescope_planes_.begin() + 1, telescope_planes_.end() - 1);
    planes_under_test_ = dia_planes_;
  } else if(at_step_ == 3){
    tel::print_banner("PART IV: Align Sil DUT (w. tracking)");
    if (sil_dut_roc_ != -1) {
      planes_to_align_.push_back(uint16_t(sil_dut_roc_));
      planes_under_test_ = vector<uint16_t>(dia_planes_.begin(), dia_planes_.end());
    }
  } else if(at_step_ == 4) {
    tel::print_banner("PART V: Align DUTs");
    planes_to_align_ = vector<uint16_t>(dia_planes_.begin(), dia_planes_.end());
    planes_under_test_ = dia_planes_;
  }
  cout << "Planes to align: " << tel::to_string(planes_to_align_) << endl;
  cout << "Planes under test : " << tel::to_string(planes_under_test_) << endl;
}

void Alignment::SetNextAlignmentStep() {

  alignment_finished_ = (at_step_ == 4) or (at_step_ == 2 and align_only_telescope_);
  at_step_++;
}

void Alignment::EventLoop(const std::vector<uint16_t> & planes) {
  ProgressBar->reset();
  ResetHistograms(); /** Reset residual histograms */
  for (uint32_t i_event = 0; FR->GetNextEvent() >= 0; ++i_event) {
    if (i_event >= max_event_number_) { break; }
    ++*ProgressBar; /** print progress */

    if (not FR->HaveOneCluster(telescope_planes_)) { continue; } /** proceed only if all telescope planes have one track */
    PLTTrack * Track = FR->Track(0);

    for (auto i_plane: planes) {
      PLTPlane * Plane = FR->Plane(i_plane);
      if (Plane->NClusters() != 1) { continue; } /** proceed only, if the plane has one cluster */

      PLTCluster * Cluster = Plane->Cluster(0);
      /** if the plane is under test, calculate the residuals otherwise just take them from the track */
      pair<float, float> dR = (FR->IsPlaneUnderTest(i_plane)) ? Track->GetResiduals(*Cluster, *FR->GetAlignment()) : Track->GetResiduals(i_plane);
      const float res_thresh = .5; // 5mm
      if (sqrt(pow(dR.first, 2) + pow(dR.second, 2)) >= res_thresh) { continue; } /** only proceed if the residual is smaller than 5mm */

      /** the point (dX, dY) has to be within the ellipsis of a=n*sig_x and b=n*sig_y: x²/a² + y²/b² <= 1 */
      if (fdX.at(i_plane).second > 0 and fdY.at(i_plane).second > 0) {  // only check if the mean residual is not 0
        if (sqrt(pow(dR.first / fdX.at(i_plane).second, 2) + pow(dR.second / fdY.at(i_plane).second, 2)) > n_sigma_) { continue; } }

      hResidual[i_plane].Fill(dR.first, dR.second); // dX vs dY
      hResidualXdY[i_plane].Fill(Cluster->LX(), dR.second); // X vs dY
      hResidualYdX[i_plane].Fill(Cluster->LY(), dR.first); // Y vs dX
    }
  }
  for (auto i_plane: planes) {
    fdX.at(i_plane) = make_pair(hResidual.at(i_plane).GetMean(1), hResidual.at(i_plane).GetRMS(1));
    fdY.at(i_plane) = make_pair(hResidual.at(i_plane).GetMean(2), hResidual.at(i_plane).GetRMS(2));

    pair<pair<float, float>, pair<float, float>> range = GetFitRange(i_plane);
    TF1 fX = TF1("fX", "pol1"), fY = TF1("fY", "pol1");
    hResidualXdY[i_plane].Fit(&fX, "q", "", range.first.first, range.first.second);
    hResidualYdX[i_plane].Fit(&fY, "q", "", range.second.first, range.second.second);
    fdA.at(i_plane) = make_pair(atan(fX.GetParameter(1)), atan(fY.GetParameter(1)));
  }
} // end EventLoop

int Alignment::Align() {

  InitGraphs();

  for (int i_align(0); i_align < maximum_steps_; i_align++) {
    cout << "BEGIN ITERATION " << i_align + 1 << " OUT OF " << maximum_steps_ << endl;
    now_ = clock();

    FR->ResetFile();
    FR->SetAllPlanes();
    if (not planes_under_test_.empty()) { FR->SetPlanesUnderTest(planes_under_test_); }

    n_sigma_ = ReduceSigma(i_align);

    EventLoop(ordered_planes_); /** Loop over all events and fill histograms */
    cout << Form("\nLoop duration: %2.1f", (clock() - now_) / CLOCKS_PER_SEC) << endl;

    for (auto roc: planes_to_align_) {
      FR->GetAlignment()->AddToLR(1, roc, (fdA.at(roc).first - fdA.at(roc).second) / 2); // take average of XdY and YdX
      FR->GetAlignment()->AddToLX(1, roc, fdX.at(roc).first);
      FR->GetAlignment()->AddToLY(1, roc, fdY.at(roc).first);
      const float cm2um = 1e4;
      g_res_mean_.at(roc)->SetPoint(i_align, i_align, cm2um * sqrt(pow(fdX.at(roc).first, 2) + pow(fdY.at(roc).first, 2)));
      g_res_angle_.at(roc)->SetPoint(i_align, i_align, fabs((fdA.at(roc).first - fdA.at(roc).second) / 2));
    }
    for (auto ipl:ordered_planes_) { SaveHistograms(ipl, i_align, at_step_); }

    CalcMaxResiduals();
    PrintResiduals(planes_to_align_);
    cout << "END ITERATION " << i_align + 1 << " OUT OF " << maximum_steps_ << endl << endl;

    if (at_step_ == 1 and i_align == 0) { ResizeHistograms(); }  // resize histograms after first real tracking

    /** Stopping criteria: max_res/angle < res/angle_thresh or the change wrt to the previous step is smaller than delta_fac * thresh*/
    if (at_step_ == 0) { // step 0 is only a translation without fit...
      break;
    } else if (GetMaxRes() < res_thresh_ and GetMaxAngle() < angle_thresh_) {
      tel::warning(Form("STOPPING ALIGNMENT:  max residual < %1.1e and max angle < %1.1e\n", res_thresh_, angle_thresh_));
      break;
    } else if (delta_max_res_.first < delta_fac_ * res_thresh_ and delta_max_res_.second < delta_fac_ * angle_thresh_) {
      tel::warning(Form("STOPPING ALIGNMENT: res correction \u0394 < %1.1e and angle correction \u0394 < %1.1e\n", delta_fac_ * res_thresh_, delta_fac_ * angle_thresh_));
      break;
    }
  } // end alignment loop

  SaveAllHistograms();
  SaveGraphs();

  PrintAlignment();
  FR->GetAlignment()->WriteAlignmentFile(telescope_id_, n_planes_);
  return 0;
}

float Alignment::ReduceSigma(uint16_t step) const {
  /** reduce n_sigma_ with each iteration until only an ellipse of "3sigma" is used to exclude residual outliers.
   *  the function n(x) = a + b / (1 + exp((x-max_steps/7) * (25 / max_steps))) is tailored to fall from b to a within ~25% of the max_steps */
  const float off_par(7), stretch_par(25); // tailored values
  return min_sigma_ + max_sigma_ / (1. + exp((step - (maximum_steps_ / off_par)) * (stretch_par / maximum_steps_)));
}

void Alignment::PrintAlignment() {

  cout << " TEL  CH ROC  RZ                             X               Y               Z" << endl;
  cout << FR->GetAlignment()->GetAlignment(telescope_id_, n_planes_) << endl;
}

void Alignment::PrintResiduals(const vector<uint16_t> & planes) {

  stringstream ss;
  ss << "RESIDUAL:  X         Y         aX        aY\n";
  for (auto i_plane: planes) {
    ss << " PLANE " << i_plane << ": " << Form("%+1.6f %+1.6f %+1.6f %+1.6f", fdX.at(i_plane).first, fdY.at(i_plane).first, fdA.at(i_plane).first,
                                               fdA.at(i_plane).second) << "\n";
  }
  ss << endl;
  ss << Form("Max residual:    %4.2e (stop at %2.1e)", GetMaxRes(), res_thresh_) << endl;
  ss << Form("Max angle:       %4.2e (stop at %2.1e)", GetMaxAngle(), angle_thresh_);
  if (last_max_res_ != delta_max_res_) {
    ss << endl;
    ss << Form("Max residual \u0394:  %4.2e (stop at %2.1e)", delta_max_res_.first, delta_fac_ * res_thresh_) << endl;
    ss << Form("Max angle \u0394:     %4.2e (stop at %2.1e)", delta_max_res_.second, delta_fac_ * angle_thresh_);
  }
  tel::print_banner(ss.str(), '-');
}

std::vector<uint16_t> Alignment::GetOrderedPlanes() {
  /** @returns the planes in the correct order wrt the beam, using the z-positions from the alignment file */
  vector<float> z_pos;
  for (unsigned i_plane(0); i_plane < n_planes_; i_plane++) {
    z_pos.emplace_back(FR->GetAlignment()->LZ(1, i_plane));
  }

  vector<uint16_t> tmp;
  for (unsigned i_plane(0); i_plane < n_planes_; i_plane++){
    auto index = distance(z_pos.begin(), min_element(z_pos.begin(), z_pos.end()));
    tmp.emplace_back(index);
    z_pos.at(index) = *max_element(z_pos.begin(), z_pos.end()) + 1;  // so it will be the highest now
  }
  cout << "Ordered planes: " << tel::to_string(tmp) << endl;
  return tmp;
}

std::vector<uint16_t> Alignment::GetTelescopePlanes(){
  /** @returns only the telescope planes (not DUT) */
  vector<uint16_t> tmp = vector<uint16_t>(ordered_planes_.begin(), ordered_planes_.begin() + n_telescope_planes / 2);
  tmp.insert(tmp.end(), ordered_planes_.end() - n_telescope_planes / 2, ordered_planes_.end());
  cout << "Telescope planes: " << tel::to_string(tmp) << endl;
  return tmp;
}

std::vector<uint16_t> Alignment::GetDiamondPlanes() {
  /** @returns the diamond pixel planes */
  uint16_t n_dia_planes = ordered_planes_.size() - n_telescope_planes;
  vector<uint16_t> tmp = {};
  if (n_dia_planes > 0) {
    tmp = vector<uint16_t>(ordered_planes_.begin() + n_telescope_planes / 2, ordered_planes_.begin() + n_telescope_planes / 2 + n_dia_planes);
    if (sil_dut_roc_ != -1){
      tmp.erase(find(tmp.begin(), tmp.end(), sil_dut_roc_));
    }
  }
  cout << "Diamond planes: " << (tmp.empty() ? "None" : tel::to_string(tmp)) << endl;
  return tmp;
}

void Alignment::InitHistograms() {
  for (uint8_t i_plane(0); i_plane != n_planes_; ++i_plane){
    hResidual.emplace_back(TH2F(Form("Residual_ROC%i", i_plane),    Form("Residual_ROC%i",    i_plane), 200, -1, 1, 300, -1, 1));
    hResidualXdY.emplace_back(TProfile(Form("ResidualXdY_ROC%i", i_plane), Form("ResidualXdY_ROC%i", i_plane), 135, -.5, .5, -1, 1));
    hResidualYdX.emplace_back(TProfile(Form("ResidualYdX_ROC%i", i_plane), Form("ResidualYdX_ROC%i", i_plane), 201, -.5, .5, -1, 1));
  }
}

void Alignment::ResizeHistograms() {

  float xmax = fdX.at(inner_planes_.back()).second * (max_sigma_ + 1);
  float ymax = fdY.at(inner_planes_.back()).second * (max_sigma_ + 1);
  for (auto ipl: telescope_planes_) {
    hResidual.at(ipl) = TH2F(Form("h%i", ipl), Form("Residual_ROC%i", ipl), 400, -xmax, xmax, 600, -ymax, ymax);
  }
  for (auto ipl: dia_planes_) {
    hResidual.at(ipl) = TH2F(Form("h%i", ipl), Form("Residual_ROC%i", ipl), 200, -xmax * 2, xmax * 2, 300, -ymax * 2, ymax * 2);
  }
}

void Alignment::ResetHistograms() {
  for (uint8_t i_plane(0); i_plane != n_planes_; ++i_plane){
    hResidual.at(i_plane).Reset();
    hResidualXdY.at(i_plane).Reset();
    hResidualYdX.at(i_plane).Reset();
  }
}

void Alignment::SaveHistograms(unsigned i_plane, int ind, int16_t step) {
  TString sub_dir = (ind != -1) ? Form("ResROC%i/", i_plane) : "";
  TString step_dir = (step != -1) ? Form("step%i/", step) : "";
  string suffix = (ind != -1) ? Form("_%02i", ind) : "";
  TString save_dir = plots_dir_ + "/" + step_dir + sub_dir;
  gSystem->mkdir(save_dir, true);
  TCanvas Can;
  Can.cd();
  // 2D Residuals
  hResidual[i_plane].SetContour(1024);
  Can.SetGridx(); Can.SetGridy();
  hResidual[i_plane].Draw("colz");
  FormatHistogram(&hResidual[i_plane], "dX [cm]", 1, "dY [cm]", 1.3, -1, 1, -1, 1);
  Can.SaveAs(save_dir + Form("%s%s%s", hResidual[i_plane].GetTitle(), suffix.c_str(), file_type_.c_str()));
  // Residual X-Projection
  gStyle->SetOptStat(1111);
  auto p_x = hResidual[i_plane].ProjectionX();
  FormatHistogram(p_x, "dX [cm]", 1, "Number of Entries", 1.3);
  p_x->Draw();
  Can.SaveAs(save_dir + Form("%s_X%s%s", hResidual[i_plane].GetTitle(), suffix.c_str(), file_type_.c_str()));
  // Residual Y-Projection
  auto p_y = hResidual[i_plane].ProjectionY();
  FormatHistogram(p_y, "dY [cm]", 1, "Number of Entries", 1.3);
  p_y->Draw();
  Can.SaveAs(save_dir + Form("%s_Y%s%s", hResidual[i_plane].GetTitle(), suffix.c_str(), file_type_.c_str()));
  // 2D Residuals X/dY
//  hResidualXdY[i_plane].SetContour(1024);
  hResidualXdY[i_plane].GetXaxis()->SetRangeUser(-0.5, 0.5);
  hResidualXdY[i_plane].GetYaxis()->SetRangeUser(-0.1, 0.1);
  hResidualXdY[i_plane].Draw();
  Can.SetGridx(); Can.SetGridy();
  Can.SaveAs(save_dir + Form("%s%s%s", hResidualXdY[i_plane].GetTitle(), suffix.c_str(), file_type_.c_str()));
  // 2D Residuals Y/dX
//  hResidualYdX[i_plane].SetContour(1024);
  hResidualYdX[i_plane].GetXaxis()->SetRangeUser(-0.5, 0.5);
  hResidualYdX[i_plane].GetYaxis()->SetRangeUser(-0.1, 0.1);
  hResidualYdX[i_plane].Draw();
  Can.SetGridx(); Can.SetGridy();
  Can.SaveAs(save_dir + Form("%s%s%s", hResidualYdX[i_plane].GetTitle(), suffix.c_str(), file_type_.c_str()));
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

  for (auto i_plane: planes_to_align_){
    SaveHistograms(i_plane, ind);
  }
}

float Alignment::CalcRes(uint16_t roc_n) {
  /** return total mean residual in the x-y-plane */
  return sqrt(pow(fdX.at(roc_n).first, 2) + pow(fdY.at(roc_n).first, 2));
}

void Alignment::CalcMaxResiduals() {
  /** calculates the max angles and residuals and pushes them in the vectors */
  float res(GetMaxRes()), angle(GetMaxAngle());
  delta_max_res_ = make_pair(fabs(last_max_res_.first - res), fabs(last_max_res_.second - angle));
  last_max_res_ = make_pair(res, angle);
}

float Alignment::GetMaxAngle() {
  /** :returns: the angle in x of the align_plane with largest angle */
  vector<float> tmp;
  for (auto i_pl: planes_to_align_) { tmp.push_back(fabs(fdA.at(i_pl).first)) ; }
  return *max_element(tmp.begin(), tmp.end());
}

float Alignment::GetMaxRes() {
  /** :returns: the residual in x-y of the align_plane with the largest residual */
  vector<float> tmp;
  for (auto i_pl: planes_to_align_) { tmp.push_back(CalcRes(i_pl)) ; }
  return *max_element(tmp.begin(), tmp.end());
}

void Alignment::ConfigROOT() {

  gROOT->ProcessLine("gErrorIgnoreLevel = kError;");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(53);
}

void Alignment::InitGraphs() {

  g_res_mean_.clear();
  g_res_angle_.clear();

  for(uint16_t i_pl = 0; i_pl < n_planes_; i_pl++){
    g_res_mean_.emplace_back(new TGraph());
    g_res_mean_[i_pl]->SetNameTitle(Form("gr%d", i_pl), Form("ResMean_Roc%d", i_pl));
    g_res_angle_.emplace_back(new TGraph());
    g_res_angle_[i_pl]->SetNameTitle(Form("ga%d", i_pl), Form("ResAngle_Roc%d", i_pl));
  }
}

void Alignment::SaveGraphs() {

  for(unsigned roc = 0; roc < n_planes_; roc++){
    gSystem->mkdir(plots_dir_, true);
    TCanvas Can;
    Can.cd();
    Can.SetGridx();
    Can.SetTickx();
    Can.SetGridy();
    Can.SetTicky();
    Can.SetLogy();
    g_res_mean_.at(roc)->GetXaxis()->SetTitle("Iteration");
    g_res_mean_.at(roc)->GetYaxis()->SetTitle("R [#mum]");
    g_res_mean_.at(roc)->GetYaxis()->SetTitleOffset(1.5);
    g_res_mean_.at(roc)->Draw("AL");
    TString fileNameCan = plots_dir_ + "/" + + Form("step%i", at_step_) + "/" + g_res_mean_.at(roc)->GetTitle();
    Can.SaveAs(Form("%s.root", fileNameCan.Data()));
    Can.SaveAs(Form("%s.png", fileNameCan.Data()));
    TCanvas Cana;
    Cana.cd();
    Cana.SetGridx();
    Cana.SetTickx();
    Cana.SetGridy();
    Cana.SetTicky();
    Cana.SetLogy();
    g_res_angle_.at(roc)->GetXaxis()->SetTitle("Iteration");
    g_res_angle_.at(roc)->GetYaxis()->SetTitle("dA [rad]");
    g_res_angle_.at(roc)->GetYaxis()->SetTitleOffset(1.5);
    g_res_angle_.at(roc)->Draw("AL");
    TString fileNameCana = plots_dir_ + "/" + Form("step%i", at_step_) + "/" + g_res_angle_.at(roc)->GetTitle();
    Cana.SaveAs(Form("%s.root", fileNameCana.Data()));
    Cana.SaveAs(Form("%s.png", fileNameCana.Data()));
  }
  cout << "\nSaved plots to: " << plots_dir_ << endl;
}

pair<pair<float, float>, pair<float, float>> Alignment::GetFitRange(uint16_t plane) {
  /** set range of profiles so that there are at least 50 entries in every bin */
  const float min = .05;
  vector<float> values;
  for (const auto & p: {hResidualXdY.at(plane), hResidualYdX.at(plane)}) {
    float xmin(p.GetBinContent(1)), xmax(p.GetBinContent(p.GetNbinsX()));
    for (auto ibin(1); ibin < p.GetNbinsX(); ibin++) { if (p.GetBinEntries(ibin) > min * p.GetEntries()) { xmin = p.GetBinCenter(ibin); break; } }
    for (auto ibin(p.GetNbinsX()); ibin > 1; ibin--) { if (p.GetBinEntries(ibin) > min * p.GetEntries()) { xmax = p.GetBinCenter(ibin); break; } }
    values.emplace_back(xmin);
    values.emplace_back(xmax);
  }
  return make_pair(make_pair(values.at(0), values.at(1)), make_pair(values.at(0), values.at(1)));
}

