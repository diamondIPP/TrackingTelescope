#include "PSIFileReader.h"

#include <iostream>
#include <string>
#include <cstdint>
#include <cstdlib>

#include "TGraph.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLine.h"
#include "GetNames.h"

using namespace std;

/** ============================
 CONSTRUCTOR
 =================================*/
PSIFileReader::PSIFileReader(bool track_only_telescope):
  PLTTracking(GetNPlanes(), track_only_telescope),
  fGainCal(fNPlanes, UseExternalCalibrationFunction()) {

    /** Set and read in gain calibration files */
    for (int i_roc=0; i_roc != fNPlanes; i_roc++) {
      string file_name = GetCalibrationPath() + Form("ROC%i.txt", i_roc);
      fCalibrationFile.emplace_back(file_name);
      fRawCalibrationFile.emplace_back(GetCalibrationPath() + Form("ph_Calibration_C%i.dat", i_roc));
      fGainCal.ReadGainCalFile(file_name, i_roc);
    }
    /** read-in additional files if we want to use the GainInterpolator */
    if (UseGainInterpolator()) {
      for (int i_roc=0; i_roc != fNPlanes; i_roc++) {
        fGainInterpolator.ReadFile(fRawCalibrationFile[i_roc], i_roc);
      }
    }

  fAlignment.ReadAlignmentFile(GetAlignmentFilename());// TODO: DA: make condition to skip this in case there is not analysis but only alignment?
  SetTrackingAlignment(&fAlignment);

  if(!trackOnlyTelescope) {
      std::cout << "Setting reader's tracking algorithm to: All plane " << std::endl;
      SetTrackingAlgorithm(PLTTracking::kTrackingAlgorithm_6PlanesHit);
  } else {
    std::cout <<"Setting reader's tracking algorithm to: ETH" << std::endl;
    SetTrackingAlgorithm(PLTTracking::kTrackingAlgorithm_ETH);
  }
}


size_t PSIFileReader::NHits ()
{
  return fHits.size();
}


PLTHit* PSIFileReader::Hit (size_t const i)
{
  if (i < fHits.size()) {
    return fHits[i];
  }
  std::cerr << "ERROR: PSIFileReader::Hit asking for a hit outside of range." << std::endl;
  throw;
}


void PSIFileReader::Clear()
{
  for (size_t i = 0; i != fHits.size(); ++i) {
    if (fHits[i]) {
      delete fHits[i];
    }
  }
  fHits.clear();
  fPlaneMap.clear();
  fPlanes.clear();
  fTracks.clear();

  return;
}


void PSIFileReader::AddToPixelMask(int ch, int roc, int col, int row)
{
  fPixelMask.insert( ch*100000 + roc*10000 + col*100 + row );
}


void PSIFileReader::ReadPixelMask (std::string const InFileName)
{
  std::cout << "PLTBinaryFileReader::ReadPixelMask reading file: " << InFileName << std::endl;

  std::ifstream InFile(InFileName.c_str());
  if (!InFile.is_open()) {
    std::cerr << "ERROR: cannot open PixelMask file: " << InFileName << std::endl;
    throw;
  }

  // Loop over header lines in the input data file
  for (std::string line; std::getline(InFile, line); ) {
    if (line == "") {
      break;
    }
  }

  for (std::string line; std::getline(InFile, line); ) {
    int ch=0, roc=0, col=0, row=0;
    std::istringstream linestream;
    linestream.str(line);

    if (line.find("col") != std::string::npos){
      linestream >> ch >> roc >> col;
      for (uint8_t i_row = 0; i_row < PLTU::NROW; i_row++)
        fPixelMask.insert(ch * 100000 + roc * 10000 + col * 100 + i_row);
    }
    else if (line.find("row") != std::string::npos){
      linestream >> ch >> roc >> row;
      for (uint8_t i_col = 0; i_col < PLTU::NCOL; i_col++)
        fPixelMask.insert(ch * 100000 + roc * 10000 + i_col * 100 + row);
    }
    else {
      linestream >> ch >> roc >> col >> row;
      fPixelMask.insert( ch*100000 + roc*10000 + col*100 + row );
    }
  }

  return;
}


bool PSIFileReader::IsPixelMasked (int const ChannelPixel)
{
  if (fPixelMask.count(ChannelPixel)) {
    return true;
  }
  return false;
}


void PSIFileReader::DrawTracksAndHits (std::string const Name)
{
    size_t const NC = NClusters();
    size_t const NT = NTracks();

    vector<float> CX(NC);
    vector<float> CY(NC);
    vector<float> CZ(NC);

  // Workaround for:
  //  TLine Line[3][NT];
  // which does not work in CLANG
    std::vector< std::vector< TLine > > Line;
    for (int i=0; i<3; i++){
        std::vector <TLine> tmp(NT);
        Line.push_back( tmp );
    }

    vector<TH2F*> HistCharge(NPlanes());
//    TH2F * HistCharge[NPlanes()];
    for (uint8_t i = 0; i != NPlanes() ; ++i) {
        TString Name;
        Name.Form("ChargeMap_ROC%i", i);
        HistCharge[i] = new TH2F(Name, Name, PLTU::NCOL, PLTU::FIRSTCOL, PLTU::LASTCOL+1, PLTU::NROW, PLTU::FIRSTROW, PLTU::LASTROW);
        HistCharge[i]->GetZaxis()->SetRangeUser(0, 50000);
    }

    vector<TH2F*> HistChargeUnclustered(NPlanes());
//    TH2F * HistChargeUnclustered[NPlanes()];
    for (uint8_t i = 0; i != NPlanes(); ++i) {
        TString Name;
        Name.Form("ChargeMapUnclustered_ROC%i", i);
        HistChargeUnclustered[i] = new TH2F(Name, Name, PLTU::NCOL, PLTU::FIRSTCOL, PLTU::LASTCOL+1, PLTU::NROW, PLTU::FIRSTROW, PLTU::LASTROW);
        HistChargeUnclustered[i]->GetZaxis()->SetRangeUser(0, 50000);
    }

    int j = 0;
    for (size_t ip = 0; ip != NPlanes(); ++ip) {
        PLTPlane* P = Plane(ip);
        for (size_t ih = 0; ih != P->NHits(); ++ih) {
            PLTHit* H = P->Hit(ih);
            ++j;

            HistCharge[ip]->SetBinContent(H->Column() + 1 - PLTU::FIRSTCOL, H->Row() + 1 - PLTU::FIRSTROW, H->Charge());
        }
        for (size_t ih = 0; ih != P->NUnclusteredHits(); ++ih) {
            PLTHit* H = P->UnclusteredHit(ih);
            HistChargeUnclustered[ip]->SetBinContent(H->Column() + 1 - PLTU::FIRSTCOL, H->Row() + 1 - PLTU::FIRSTROW, H->Charge());
        }
    }

    float zMax(0), zMin(20);
    int jc = 0;
    for (size_t ip = 0; ip != NPlanes(); ++ip) {
        PLTPlane* P = Plane(ip);
        for (size_t ic = 0; ic != P->NClusters(); ++ic) {
            PLTCluster* C = P->Cluster(ic);
            CX[jc] = C->TX();
            CY[jc] = C->TY();
            CZ[jc] = C->TZ();
            if (C->TZ() > zMax) zMax = C->TZ();
            if (C->TZ() < zMin) zMin = C->TZ();
            ++jc;
        }
    }

    std::vector<PLTHit*> UsedHits;

    for (uint8_t i = 0; i != NT; ++i) {
        PLTTrack* T = fTracks[i];


        // XZ
        Line[0][i].SetX1(zMin);
        Line[0][i].SetX2(zMax);
        Line[0][i].SetY1(T->TX(zMin));
        Line[0][i].SetY2(T->TX(zMax));
        Line[0][i].SetLineColor(i+1);

        // YZ
        Line[1][i].SetX1(zMin);
        Line[1][i].SetX2(zMax);
        Line[1][i].SetY1(T->TY(zMin));
        Line[1][i].SetY2(T->TY(zMax));
        Line[1][i].SetLineColor(i+1);

        // XY
        Line[2][i].SetX1(T->TX(zMin));
        Line[2][i].SetX2(T->TX(zMax));
        Line[2][i].SetY1(T->TY(zMin));
        Line[2][i].SetY2(T->TY(zMax));
        Line[2][i].SetLineColor(i+1);

    //printf("XY 0 7: %9.3f %9.3f   %9.3f %9.3f\n", T->TX(0), T->TY(0), T->TX(7), T->TY(7));
  }

  //for (int i = 0; i != 3; ++i) {
  //  for (int j = 0; j != NT; ++j) {
  //    Line[i][j].SetLineColor(4);
  //  }
  //}

  TCanvas C("TelescopeTrack", "TelescopeTrack", 800, 800);;
  C.Divide(3, 4);

  C.cd(1);
  TGraph gXZ(NC, &CZ[0], &CX[0]);
  gXZ.SetTitle("");
  gXZ.GetXaxis()->SetTitle("Z (cm)");
  gXZ.GetYaxis()->SetTitle("X (cm)");
  gXZ.GetXaxis()->SetTitleSize(0.06);
  gXZ.GetYaxis()->SetTitleSize(0.08);
  gXZ.GetXaxis()->SetTitleOffset(0.7);
  gXZ.GetYaxis()->SetTitleOffset(0.5);
  gXZ.SetMarkerColor(40);
  gXZ.GetXaxis()->SetLimits(zMin - 0.5, zMax + 1);
  gXZ.SetMinimum(-0.5);
  gXZ.SetMaximum( 0.5);
  if (NC) {
    gXZ.Draw("A*");
  }
  for (size_t i = 0; i != NT; ++i) {
    Line[0][i].Draw();
  }

  C.cd(2);
  TGraph gYZ(NC, &CZ[0], &CY[0]);
  gYZ.SetTitle("");
  gYZ.GetXaxis()->SetTitle("Z (cm)");
  gYZ.GetYaxis()->SetTitle("Y (cm)");
  gYZ.GetXaxis()->SetTitleSize(0.06);
  gYZ.GetYaxis()->SetTitleSize(0.08);
  gYZ.GetXaxis()->SetTitleOffset(0.7);
  gYZ.GetYaxis()->SetTitleOffset(0.5);
  gYZ.SetMarkerColor(40);
  gYZ.GetXaxis()->SetLimits(zMin - 0.5, zMax + 1);
  gYZ.SetMinimum(-0.5);
  gYZ.SetMaximum( 0.5);
  if (NC) {
    gYZ.Draw("A*");
  }
  for (size_t i = 0; i != NT; ++i) {
    Line[1][i].Draw();
  }

  //TVirtualPad* Pad = C.cd(3);
  //Pad->DrawFrame(-30, -30, 30, 30);
  C.cd(3);
  TGraph gXY(NC, &CX[0], &CY[0]);
  gXY.SetTitle("");
  gXY.GetXaxis()->SetTitle("X (cm)");
  gXY.GetYaxis()->SetTitle("Y (cm)");
  gXY.GetXaxis()->SetTitleSize(0.06);
  gXY.GetYaxis()->SetTitleSize(0.08);
  gXY.GetXaxis()->SetTitleOffset(0.7);
  gXY.GetYaxis()->SetTitleOffset(0.5);
  gXY.SetMarkerColor(40);
  gXY.GetXaxis()->SetLimits(-0.5, 0.5);
  gXY.SetMinimum(-0.5);
  gXY.SetMaximum( 0.5);
  if (NC) {
    gXY.Draw("A*");
  }
  for (size_t i = 0; i != NT; ++i) {
    Line[2][i].Draw();
  }

    for (uint8_t iHist = 0; iHist != NPlanes(); iHist++){
        C.cd(iHist + 4);
        HistCharge[iHist]->Draw("colz");
    }
  //C.cd(3);
  //HistChargeUnclustered[0]->Draw("colz");
  //C.cd(6);
  //HistChargeUnclustered[1]->Draw("colz");
  //C.cd(9);
  //HistChargeUnclustered[2]->Draw("colz");

  C.SaveAs(Name.c_str());

    for (uint8_t i = 0; i != NPlanes(); ++i) {
        delete HistCharge[i];
        delete HistChargeUnclustered[i];
    }

  return;
}

