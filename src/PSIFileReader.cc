#include "PSIFileReader.h"

#include <iostream>
#include <string>
#include <stdint.h>
#include <stdlib.h>

#include "TGraph.h"
#include "TString.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include "TMarker.h"
#include "TLine.h"

using namespace std;

/** ============================
 CONSTRUCTOR
 =================================*/
PSIFileReader::PSIFileReader (string const CalibrationList, string const AlignmentFileName,
    int const nrocs, bool const useGainInterpolator, bool const useExternalCalibrationFunction):
        PLTTracking(nrocs), NMAXROCS(nrocs), fGainCal(nrocs, useExternalCalibrationFunction),
        fUseGainInterpolator(useGainInterpolator){

   /** Initialize fCalibrationFile and fRawCalibrationFile with empty strings */
    for (int i_roc=0; i_roc != NMAXROCS; i_roc++) {
        fCalibrationFile.push_back("");
        fRawCalibrationFile.push_back("");
    }
    std::ifstream fCL(CalibrationList.c_str());
    if (!fCL.is_open()) {
        std::cerr << "ERROR: cannot open calibration list file: " << CalibrationList << std::endl;
        throw;
    }
    fCL >> fBaseCalDir;
    for (int i_roc=0; i_roc != NMAXROCS; i_roc++) fCL >> fCalibrationFile[i_roc];

    /** look for additional files if we want to use the GainInterpolator */
    if (useGainInterpolator) {
        for (int i_roc=0; i_roc != NMAXROCS; i_roc++)
            fCL >> fRawCalibrationFile[i_roc];
    }

    for (int i_roc=0; i_roc != NMAXROCS; i_roc++)
        fGainCal.ReadGainCalFile(fBaseCalDir + "/" + fCalibrationFile[i_roc], i_roc);

    /** read-in additional files if we want to use the GainInterpolator */
    if (useGainInterpolator) {
        for (int i_roc=0; i_roc != NMAXROCS; i_roc++)
            fGainInterpolator.ReadFile(fBaseCalDir + "/" + fRawCalibrationFile[i_roc], i_roc);
    }

    fAlignment.ReadAlignmentFile(AlignmentFileName);
    SetTrackingAlignment(&fAlignment);

    SetTrackingAlgorithm(PLTTracking::kTrackingAlgorithm_6PlanesHit);
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
    linestream >> ch >> roc >> col >> row;

    //std::cout << "Masking: " << ch << " " << roc << " " << col << " " << row << std::endl;
    fPixelMask.insert( ch*100000 + roc*10000 + col*100 + row );
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
  int const NH = NHits();
  int const NC = NClusters();
  int const NT = NTracks();

  float X[NH];
  float Y[NH];
  float Z[NH];

  float CX[NC];
  float CY[NC];
  float CZ[NC];

  // Workaround for:
  //  TLine Line[3][NT];
  // which does not work in CLANG
  std::vector< std::vector< TLine > > Line;
  for (int i=0; i<3; i++){
    std::vector <TLine> tmp(NT);
    Line.push_back( tmp );
  }

  TH2F* HistCharge[6];
  for (int i = 0; i != 6; ++i) {
    TString Name;
    Name.Form("ChargeMap_ROC%i", i);
    HistCharge[i] = new TH2F(Name, Name, PLTU::NCOL, PLTU::FIRSTCOL, PLTU::LASTCOL+1, PLTU::NROW, PLTU::FIRSTROW, PLTU::LASTROW);
    HistCharge[i]->GetZaxis()->SetRangeUser(0, 50000);
  }

  TH2F* HistChargeUnclustered[6];
  for (int i = 0; i != 6; ++i) {
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
      X[j] = H->TX();
      Y[j] = H->TY();
      Z[j] = H->TZ();
      ++j;

      HistCharge[ip]->SetBinContent(H->Column() + 1 - PLTU::FIRSTCOL, H->Row() + 1 - PLTU::FIRSTROW, H->Charge());
    }
    for (size_t ih = 0; ih != P->NUnclusteredHits(); ++ih) {
      PLTHit* H = P->UnclusteredHit(ih);
      HistChargeUnclustered[ip]->SetBinContent(H->Column() + 1 - PLTU::FIRSTCOL, H->Row() + 1 - PLTU::FIRSTROW, H->Charge());
    }



  }
  int jc = 0;
  for (size_t ip = 0; ip != NPlanes(); ++ip) {
    PLTPlane* P = Plane(ip);
    for (size_t ic = 0; ic != P->NClusters(); ++ic) {
      PLTCluster* C = P->Cluster(ic);
      CX[jc] = C->TX();
      CY[jc] = C->TY();
      CZ[jc] = C->TZ();
      ++jc;
    }
  }

  std::vector<PLTHit*> UsedHits;

  for (int i = 0; i != NT; ++i) {
    PLTTrack* T = fTracks[i];

    // XZ
    Line[0][i].SetX1(0);
    Line[0][i].SetX2(10.16);
    Line[0][i].SetY1(T->TX(0));
    Line[0][i].SetY2(T->TX(10.16));
    Line[0][i].SetLineColor(i+1);

    // YZ
    Line[1][i].SetX1(0);
    Line[1][i].SetX2(10.16);
    Line[1][i].SetY1(T->TY(0));
    Line[1][i].SetY2(T->TY(10.16));
    Line[1][i].SetLineColor(i+1);

    // XY
    Line[2][i].SetX1(T->TX(0));
    Line[2][i].SetX2(T->TX(10.16));
    Line[2][i].SetY1(T->TY(0));
    Line[2][i].SetY2(T->TY(10.16));
    Line[2][i].SetLineColor(i+1);

    //printf("XY 0 7: %9.3f %9.3f   %9.3f %9.3f\n", T->TX(0), T->TY(0), T->TX(7), T->TY(7));
  }

  //for (int i = 0; i != 3; ++i) {
  //  for (int j = 0; j != NT; ++j) {
  //    Line[i][j].SetLineColor(4);
  //  }
  //}

  TCanvas C("TelescopeTrack", "TelescopeTrack", 800, 800);;
  C.Divide(3, 3);

  C.cd(1);
  TGraph gXZ(NC, CZ, CX);
  gXZ.SetTitle("");
  gXZ.GetXaxis()->SetTitle("Z (cm)");
  gXZ.GetYaxis()->SetTitle("X (cm)");
  gXZ.GetXaxis()->SetTitleSize(0.06);
  gXZ.GetYaxis()->SetTitleSize(0.08);
  gXZ.GetXaxis()->SetTitleOffset(0.7);
  gXZ.GetYaxis()->SetTitleOffset(0.5);
  gXZ.SetMarkerColor(40);
  gXZ.GetXaxis()->SetLimits(-0.5, 11.0);
  gXZ.SetMinimum(-0.5);
  gXZ.SetMaximum( 0.5);
  if (NC) {
    gXZ.Draw("A*");
  }
  for (int i = 0; i != NT; ++i) {
    Line[0][i].Draw();
  }

  C.cd(4);
  TGraph gYZ(NC, CZ, CY);
  gYZ.SetTitle("");
  gYZ.GetXaxis()->SetTitle("Z (cm)");
  gYZ.GetYaxis()->SetTitle("Y (cm)");
  gYZ.GetXaxis()->SetTitleSize(0.06);
  gYZ.GetYaxis()->SetTitleSize(0.08);
  gYZ.GetXaxis()->SetTitleOffset(0.7);
  gYZ.GetYaxis()->SetTitleOffset(0.5);
  gYZ.SetMarkerColor(40);
  gYZ.GetXaxis()->SetLimits(-0.5, 11.0);
  gYZ.SetMinimum(-0.5);
  gYZ.SetMaximum( 0.5);
  if (NC) {
    gYZ.Draw("A*");
  }
  for (int i = 0; i != NT; ++i) {
    Line[1][i].Draw();
  }

  //TVirtualPad* Pad = C.cd(3);
  //Pad->DrawFrame(-30, -30, 30, 30);
  C.cd(7);
  TGraph gXY(NC, CX, CY);
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
  for (int i = 0; i != NT; ++i) {
    Line[2][i].Draw();
  }

  C.cd(2);
  HistCharge[0]->Draw("colz");
  C.cd(3);
  HistCharge[1]->Draw("colz");
  C.cd(5);
  HistCharge[2]->Draw("colz");
  C.cd(6);
  HistCharge[3]->Draw("colz");
  C.cd(8);
  HistCharge[4]->Draw("colz");
  C.cd(9);
  HistCharge[5]->Draw("colz");

  //C.cd(3);
  //HistChargeUnclustered[0]->Draw("colz");
  //C.cd(6);
  //HistChargeUnclustered[1]->Draw("colz");
  //C.cd(9);
  //HistChargeUnclustered[2]->Draw("colz");

  C.SaveAs(Name.c_str());

  for (int i = 0; i != 6; ++i) {
    delete HistCharge[i];
    delete HistChargeUnclustered[i];
  }

  return;
}

