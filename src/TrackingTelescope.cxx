////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <Dean.Andrew.Hidas@cern.ch>
//
// Created on: Wed Apr  9 13:49:09 CEST 2014
//
////////////////////////////////////////////////////////////////////

#include <iostream>
#include <utility>
#include <cmath>
// numeric is needed to compile with make in certain systems
#include <numeric>
#include <algorithm>
#include <zconf.h>

#include "TLegend.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH3F.h"
#include "TProfile2D.h"
#include "TParameter.h"
#include "TInterpreter.h"

#include "PLTAnalysis.h"

#include "PSIBinaryFileReader.h"
#include "PSIRootFileReader.h"
#include "DoAlignment.h"
#include "FindPlaneErrors.h"
#include "Utils.h"
#include "GetNames.h"

#define DEBUG false

using namespace std;

template<typename T>
void FillIth(TH3F *h, int px, int py, std::vector<T> values, int i, bool fill_when_too_small = true){
/*
    Fill the i-th value of the vector into the histogram.
    px and py are the position to fill them at.

    If i is negative: count from the back, Python style!

    fill_when_too_small decides if we should fill with zero when not enough values are available,
*/
    // Count from the front
    // [0,1,2,3,...]
    if (i>=0){
        if (values.size() >= unsigned(i + 1))
            h->Fill(px, py, values[i]);
        else
            if (fill_when_too_small)
                h->Fill(px, py, 0);
    }
    // Count from the back
    // [...,-3, -2, -1]
    else{
      int abs_i = abs(i);

      if (values.size() >= unsigned(abs_i))
        h->Fill(px, py, values[values.size()-abs_i]);
      else
        if (fill_when_too_small)
            h->Fill(px, py, 0);
    }

}

// Helper function for sorting a <int, float> pair according to the float
bool PairSort( std::pair<int, float> i, std::pair<int, float> j) { return (i.second < j.second); }

std::pair<int, float> FindILowestIndexAndValue( std::vector<float> values, uint32_t i=0){

    std::vector<std::pair<int, float> > index_and_values;

    if (DEBUG){
        std::cout << "FindILowestIndexAndValue, before sort: ";
        for (uint32_t iv = 0; iv != values.size(); iv++)
            std::cout << values[iv] << " ";
        std::cout << endl;
    }

    // SORT
    for (uint32_t iv = 0; iv != values.size(); iv++)
        index_and_values.push_back(std::make_pair(iv, values[iv]));


    std::sort(index_and_values.begin(), index_and_values.end(), PairSort);

   if (DEBUG){
        std::cout << "FindILowestIndexAndValue, after sort: ";
        for (uint32_t iv = 0; iv != index_and_values.size(); iv++)
            std::cout << index_and_values[iv].second << " ";
        std::cout << endl;
    }


    if ((i+1) <= index_and_values.size()){
        if (DEBUG)
            std::cout << "FindILowestIndexAndValue, i=" << i << " " << index_and_values[i].first << " : " << index_and_values[i].second << std::endl;
        return index_and_values[i];
    }
    else
        return std::make_pair(-1, TMath::QuietNaN());

}



bool CheckEllipse(float dx, float dy, float max_dx, float max_dy){
  if ( (dx*dx/(max_dx*max_dx)) + (dy*dy/(max_dy*max_dy)) <= 1.01)
    return true;
  else
    return false;
}

void Write2DCharge( TH3* h, TCanvas * Can, float maxz, TString OutDir){
  TProfile2D * ph = h->Project3DProfile("yx");
  ph->SetAxisRange(12,38,"X");
  ph->SetAxisRange(39,80,"Y");
  ph->SetMinimum(0);
  ph->SetMaximum(maxz);
  ph->Draw("COLZ");
  ph->Write();
  Can->SaveAs( OutDir+ TString(h->GetName()) +"_profile.png");
  Can->SaveAs( OutDir+ TString(h->GetName()) +"_profile.pdf");
}


void WriteAngleHistograms( TH1 * h_before_chi2_x,
                          TH1 * h_before_chi2_y,
                          TH1 * h_after_chi2_x,
                          TH1 * h_after_chi2_y,
                          TCanvas * Can,
                          TString OutDir){

  h_before_chi2_x->SetLineColor( kRed );
  h_before_chi2_y->SetLineColor( kBlue );
  h_after_chi2_x->SetLineColor( kMagenta );
  h_after_chi2_y->SetLineColor( kBlack );

  h_before_chi2_x->SetLineStyle(1);
  h_before_chi2_y->SetLineStyle(2);
  h_after_chi2_x->SetLineStyle(3);
  h_after_chi2_y->SetLineStyle(4);

  h_before_chi2_x->SetLineWidth(2);
  h_before_chi2_y->SetLineWidth(2);
  h_after_chi2_x->SetLineWidth(2);
  h_after_chi2_y->SetLineWidth(2);

  float hmax = 1.1 * std::max( h_before_chi2_x->GetMaximum(),
                         std::max( h_before_chi2_y->GetMaximum(),
                             std::max( h_after_chi2_x->GetMaximum(),
                                h_after_chi2_y->GetMaximum())));


  h_before_chi2_x->SetAxisRange(0, hmax, "Y");
  h_before_chi2_y->SetAxisRange(0, hmax, "Y");
  h_after_chi2_x->SetAxisRange(0, hmax, "Y");
  h_after_chi2_y->SetAxisRange(0, hmax, "Y");

  h_before_chi2_x->GetXaxis()->SetTitle("Angle [rad]");
  h_before_chi2_x->GetYaxis()->SetTitle("Tracks");

  h_before_chi2_x->Draw();
  h_before_chi2_y->Draw("SAME");
  h_after_chi2_x->Draw("SAME");
  h_after_chi2_y->Draw("SAME");


  TLegend Leg(0.7, 0.5, 0.85, 0.88, "");
  Leg.SetFillColor(0);
  Leg.SetBorderSize(0);
  Leg.SetTextSize(0.04);
  Leg.AddEntry(h_before_chi2_x, "X Before Chi^{2}", "l");
  Leg.AddEntry(h_before_chi2_y, "Y Before Chi^{2}", "l");
  Leg.AddEntry(h_after_chi2_x, "X After Chi^{2}", "l");
  Leg.AddEntry(h_after_chi2_y, "Y After Chi^{2}", "l");

  Leg.Draw();

  Can->SaveAs( OutDir+ TString(h_before_chi2_x->GetName()) + ".png");

}


float GetMaximumExceptBin(TH1* h, int ibin){

  if (h->GetMaximumBin() == ibin){
    return h->GetMaximum(h->GetMaximum());
  }
  else{
    return h->GetMaximum();
  }

}

void Write1DCharge( std::vector<TH3*> hs, TCanvas *Can, TString OutDir){

  if (hs.size()!=4){
    std::cerr << "Write1DCharge needs exactly four histograms!" << std::endl;
    return;
  }

  TH1* h15 = hs[0]->Project3D("Z");
  TH1* h30 = hs[1]->Project3D("Z");
  TH1* h45 = hs[2]->Project3D("Z");
  TH1* h60 = hs[3]->Project3D("Z");

  float hmax = 1.1 * std::max( GetMaximumExceptBin(h15, 1),
                       std::max( GetMaximumExceptBin(h30, 1),
                         std::max( GetMaximumExceptBin(h45, 1),
                            GetMaximumExceptBin(h60, 1))));

  h15->SetAxisRange(0,hmax,"Y");
  h30->SetAxisRange(0,hmax,"Y");
  h45->SetAxisRange(0,hmax,"Y");
  h60->SetAxisRange(0,hmax,"Y");


  h15->SetLineColor(1);
  h30->SetLineColor(2);
  h45->SetLineColor(3);
  h60->SetLineColor(4);

  h15->SetLineWidth(2);
  h30->SetLineWidth(2);
  h45->SetLineWidth(2);
  h60->SetLineWidth(2);

  h15->GetXaxis()->SetTitle("Charge (Electrons)");
  h15->GetYaxis()->SetTitle("Number of Hits");

  TLegend Leg(0.7, 0.5, 0.90, 0.88, "");
  Leg.SetFillColor(0);
  Leg.SetBorderSize(0);
  Leg.SetTextSize(0.05);
  Leg.AddEntry(h15, "R=1", "l");
  Leg.AddEntry(h30, "R=2", "l");
  Leg.AddEntry(h45, "R=3", "l");
  Leg.AddEntry(h60, "R=4", "l");

  h15->Draw();
  h30->Draw("SAME");
  h45->Draw("SAME");
  h60->Draw("SAME");
  Leg.Draw();

  h15->Write();
  h30->Write();
  h45->Write();
  h60->Write();

  Can->SaveAs( OutDir+ TString(hs[0]->GetName()) +".png");
  Can->SaveAs( OutDir+ TString(hs[0]->GetName()) +".pdf");

}

void Write1DFraction(std::vector<TH1*> hs, TCanvas *Can, TString OutDir){

  if (hs.size()!=4){
    std::cerr << "Write1Dneeds exactly four histograms!" << std::endl;
    return;
  }


  float hmax = 1.1 * std::max( hs[0]->GetMaximum(),
                       std::max( hs[1]->GetMaximum(),
                         std::max( hs[2]->GetMaximum(),
                            hs[3]->GetMaximum())));

  hs[0]->SetAxisRange(0,hmax,"Y");
  hs[1]->SetAxisRange(0,hmax,"Y");
  hs[2]->SetAxisRange(0,hmax,"Y");
  hs[3]->SetAxisRange(0,hmax,"Y");


  hs[0]->SetLineColor(1);
  hs[1]->SetLineColor(2);
  hs[2]->SetLineColor(3);
  hs[3]->SetLineColor(4);

  hs[0]->SetLineWidth(2);
  hs[1]->SetLineWidth(2);
  hs[2]->SetLineWidth(2);
  hs[3]->SetLineWidth(2);

  hs[0]->GetXaxis()->SetTitle("Fraction of Hits in Cluster inside Radius");
  hs[0]->GetYaxis()->SetTitle("");

  TLegend Leg(0.7, 0.5, 0.90, 0.88, "");
  Leg.SetFillColor(0);
  Leg.SetBorderSize(0);
  Leg.SetTextSize(0.05);
  Leg.AddEntry(hs[0], "R=1", "l");
  Leg.AddEntry(hs[1], "R=2", "l");
  Leg.AddEntry(hs[2], "R=3", "l");
  Leg.AddEntry(hs[3], "R=4", "l");

  hs[0]->Draw();
  hs[1]->Draw("SAME");
  hs[2]->Draw("SAME");
  hs[3]->Draw("SAME");
  Leg.Draw();

  Can->SaveAs( OutDir+ TString(hs[0]->GetName()) +".png");
  Can->SaveAs( OutDir+ TString(hs[0]->GetName()) +".pdf");

}



void TestPlaneEfficiency (std::string const InFileName,
                          TFile * out_f,
                          TString const RunNumber,
                          int plane_under_test,
                          int n_events,
                          int telescopeID)
{
  /* TestPlaneEfficiency

  o) Consider one plane to be the plane under test
  o) Require exactly one hit in all other planes
  o) This gives one track
  o) Then check if a hit was registered in the plane under test (within a given
      radius around the expected passing of the track)

  */

  // Track/Hit matching distance [cm]
  float max_dr_x = 0.03;
  float max_dr_y = 0.02;

  gStyle->SetOptStat(0);

  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);

  TString const PlotsDir = "plots/";
  TString const OutDir = PlotsDir + RunNumber + "/";

  // Initialize Reader
  PSIFileReader * FR;

  if (IsROOTFile(InFileName)){
    FR = new PSIRootFileReader(InFileName, 0, false);
  }
  else{
    FR = new PSIBinaryFileReader(InFileName);
    ((PSIBinaryFileReader*) FR)->CalculateLevels(OutDir);
  }

  FR->GetAlignment()->SetErrors(telescopeID);
  FR->SetPlaneUnderTest(plane_under_test);

  // Apply Masking
  FR->ReadPixelMask(GetMaskingFilename());

  // Prepare Occupancy histograms
  // Telescope coordinates
  TH2F hOccupancyNum   = TH2F(Form("PlaneEfficiency_ROC%i", plane_under_test), "PlaneEfficiency",   52, 0, 52, 80, 0, 80);
  TH2F hOccupancyDenom = TH2F(Form("TracksPassing_ROC%i", plane_under_test), Form("TracksPassing_ROC%i",plane_under_test), 52, 0, 52, 80, 0, 80);

  // Also have a second set - sliced according to event number
  int n_slices = 5;
  int slice_size = n_events/n_slices;
  std::vector<TH2F> hOccupancyNum_eventSlices;
  std::vector<TH2F> hOccupancyDenom_eventSlices;
  for (int i=0; i != n_slices; i++){

    TH2F h_n = TH2F(   Form("Numerator_ROC%i_slice%i",plane_under_test,i), "",   52, 0, 52, 80, 0, 80);
    TH2F h_d = TH2F(   Form("Denominator_ROC%i_slice%i",plane_under_test,i), "",   52, 0, 52, 80, 0, 80);

    hOccupancyNum_eventSlices.push_back( h_n );
    hOccupancyDenom_eventSlices.push_back( h_d );
  }

  TH3F hSumCharge1 = TH3F( Form("SumCharge_ROC%i", plane_under_test),  "Sum Charge within 1-Pixel Ellipse", 52,0,52, 80,0,80,50,0,50000);
  TH3F hSumCharge2 = TH3F( Form("SumCharge2_ROC%i", plane_under_test), "Sum Charge within 2-Pixel Ellipse", 52,0,52, 80,0,80,50,0,50000);
  TH3F hSumCharge3 = TH3F( Form("SumCharge3_ROC%i", plane_under_test), "Sum Charge within 3-Pixel Ellipse", 52,0,52, 80,0,80,50,0,50000);
  TH3F hSumCharge4 = TH3F( Form("SumCharge4_ROC%i", plane_under_test), "Sum Charge within 4-Pixel Ellipse", 52,0,52, 80,0,80,50,0,50000);

  TH3F h1stCharge1 = TH3F( Form("1stCharge_ROC%i", plane_under_test),  "1st Charge within 1-Pixel Ellipse", 52,0,52, 80,0,80,50,0,50000);
  TH3F h1stCharge2 = TH3F( Form("1stCharge2_ROC%i", plane_under_test), "1st Charge within 2-Pixel Ellipse", 52,0,52, 80,0,80,50,0,50000);
  TH3F h1stCharge3 = TH3F( Form("1stCharge3_ROC%i", plane_under_test), "1st Charge within 3-Pixel Ellipse", 52,0,52, 80,0,80,50,0,50000);
  TH3F h1stCharge4 = TH3F( Form("1stCharge4_ROC%i", plane_under_test), "1st Charge within 4-Pixel Ellipse", 52,0,52, 80,0,80,50,0,50000);

  TH3F h1stCharge1ADC = TH3F( Form("1stCharge_ADC_ROC%i", plane_under_test),  "1st Charge within 1-Pixel Ellipse", 52,0,52, 80,0,80, 50,-700, -200);
  TH3F h1stCharge2ADC = TH3F( Form("1stCharge2_ADC_ROC%i", plane_under_test), "1st Charge within 2-Pixel Ellipse", 52,0,52, 80,0,80, 50,-700, -200);
  TH3F h1stCharge3ADC = TH3F( Form("1stCharge3_ADC_ROC%i", plane_under_test), "1st Charge within 3-Pixel Ellipse", 52,0,52, 80,0,80, 50,-700, -200);
  TH3F h1stCharge4ADC = TH3F( Form("1stCharge4_ADC_ROC%i", plane_under_test), "1st Charge within 4-Pixel Ellipse", 52,0,52, 80,0,80, 50,-700, -200);

  TH3F h2ndCharge1 = TH3F( Form("2ndCharge_ROC%i", plane_under_test),  "2nd Charge within 1-Pixel Ellipse", 52,0,52, 80,0,80,50,0,50000);
  TH3F h2ndCharge2 = TH3F( Form("2ndCharge2_ROC%i", plane_under_test), "2nd Charge within 2-Pixel Ellipse", 52,0,52, 80,0,80,50,0,50000);
  TH3F h2ndCharge3 = TH3F( Form("2ndCharge3_ROC%i", plane_under_test), "2nd Charge within 3-Pixel Ellipse", 52,0,52, 80,0,80,50,0,50000);
  TH3F h2ndCharge4 = TH3F( Form("2ndCharge4_ROC%i", plane_under_test), "2nd Charge within 4-Pixel Ellipse", 52,0,52, 80,0,80,50,0,50000);

  TH3F h2ndCharge1ADC = TH3F( Form("2ndCharge_ADC_ROC%i", plane_under_test),  "2nd Charge within 1-Pixel Ellipse", 52,0,52, 80,0,80, 50, -700, -200);
  TH3F h2ndCharge2ADC = TH3F( Form("2ndCharge2_ADC_ROC%i", plane_under_test), "2nd Charge within 2-Pixel Ellipse", 52,0,52, 80,0,80, 50, -700, -200);
  TH3F h2ndCharge3ADC = TH3F( Form("2ndCharge3_ADC_ROC%i", plane_under_test), "2nd Charge within 3-Pixel Ellipse", 52,0,52, 80,0,80, 50, -700, -200);
  TH3F h2ndCharge4ADC = TH3F( Form("2ndCharge4_ADC_ROC%i", plane_under_test), "2nd Charge within 4-Pixel Ellipse", 52,0,52, 80,0,80, 50, -700, -200);

  TH1F hFractionContainted1 = TH1F(Form("FractionContained1_ROC%i", plane_under_test), "Fraction Contained", 50, 0, 1.1);
  TH1F hFractionContainted2 = TH1F(Form("FractionContained2_ROC%i", plane_under_test), "Fraction Contained", 50, 0, 1.1);
  TH1F hFractionContainted3 = TH1F(Form("FractionContained3_ROC%i", plane_under_test), "Fraction Contained", 50, 0, 1.1);
  TH1F hFractionContainted4 = TH1F(Form("FractionContained4_ROC%i", plane_under_test), "Fraction Contained", 50, 0, 1.1);


  TH1F hCS1_1stCharge = TH1F( Form("CS1_1stCharge_ROC%i",plane_under_test),   "Cluster Size 1, 1st Charge", 100, 0, 50000);
  TH1F hCS1_SumCharge = TH1F( Form("CS1_SumCharge_ROC%i",plane_under_test),   "Cluster Size 1, Sum Charge", 100, 0, 50000);

  TH1F hCS2_1stCharge = TH1F( Form("CS2_1stCharge_ROC%i",plane_under_test),   "Cluster Size 2, 1st Charge", 100, 0, 50000);
  TH1F hCS2_2ndCharge = TH1F( Form("CS2_2ndCharge_ROC%i",plane_under_test),   "Cluster Size 2, 2nd Charge", 100, 0, 50000);
  TH1F hCS2_SumCharge = TH1F( Form("CS2_SumCharge_ROC%i",plane_under_test),   "Cluster Size 2, Sum Charge", 100, 0, 50000);

  TH2F hCS2_2D = TH2F( Form("CS2_2D_ROC%i",plane_under_test),   "Cluster Size 2, 2D", 100, 0, 50000, 100, 0, 50000);

  TH1F hCS3_1stCharge = TH1F( Form("CS3_1stCharge_ROC%i",plane_under_test),   "Cluster Size 3, 1st Charge", 100, 0, 50000);
  TH1F hCS3_2ndCharge = TH1F( Form("CS3_2ndCharge_ROC%i",plane_under_test),   "Cluster Size 3, 2nd Charge", 100, 0, 50000);
  TH1F hCS3_3rdCharge = TH1F( Form("CS3_3rdCharge_ROC%i",plane_under_test),   "Cluster Size 3, 3rd Charge", 100, 0, 50000);
  TH1F hCS3_SumCharge = TH1F( Form("CS3_SumCharge_ROC%i",plane_under_test),   "Cluster Size 3, Sum Charge", 100, 0, 50000);

  TH1F hCS4_1stCharge = TH1F( Form("CS4_1stCharge_ROC%i",plane_under_test),   "Cluster Size 4, 1st Charge", 100, 0, 50000);
  TH1F hCS4_2ndCharge = TH1F( Form("CS4_2ndCharge_ROC%i",plane_under_test),   "Cluster Size 4, 2nd Charge", 100, 0, 50000);
  TH1F hCS4_3rdCharge = TH1F( Form("CS4_3rdCharge_ROC%i",plane_under_test),   "Cluster Size 4, 3rd Charge", 100, 0, 50000);
  TH1F hCS4_4thCharge = TH1F( Form("CS4_4thCharge_ROC%i",plane_under_test),   "Cluster Size 4, 4th Charge", 100, 0, 50000);
  TH1F hCS4_SumCharge = TH1F( Form("CS4_SumCharge_ROC%i",plane_under_test),   "Cluster Size 4, Sum Charge", 100, 0, 50000);


  TH3F hClusterSize       = TH3F( Form("ClusterSize_ROC%i", plane_under_test), "Cluster Size", 52,0,52, 80,0,80,11,-0.5,10.5);

  TH1F hdtx = TH1F( Form("SinglePlaneTestDX_ROC%i",plane_under_test),   "SinglePlaneTest_DX",   100, -0.2, 0.2 );
  TH1F hdty = TH1F( Form("SinglePlaneTestDY_ROC%i",plane_under_test),   "SinglePlaneTest_DY",   100, -0.2, 0.2 );
  TH1F hdtr = TH1F( Form("SinglePlaneTestDR_ROC%i",plane_under_test),   "SinglePlaneTest_DR",   100, 0, 0.4 );

  TH1F hDrSecondCluster = TH1F(Form("DeltaRSecondCluster_ROC%i", plane_under_test), "#Delta R Second Cluster", 50, -2., 20);

  TH1F hChi2  = TH1F( Form("SinglePlaneTestChi2_ROC%i",plane_under_test),   "SinglePlaneTest_Chi2",    200, 0, 50 );
  TH1F hChi2X = TH1F( Form("SinglePlaneTestChi2X_ROC%i",plane_under_test),  "SinglePlaneTest_Chi2X",   100, 0, 20 );
  TH1F hChi2Y = TH1F( Form("SinglePlaneTestChi2Y_ROC%i",plane_under_test),  "SinglePlaneTest_Chi2Y",   100, 0, 20 );


  hChi2X.GetXaxis()->SetTitleSize(0.06);
  hChi2X.GetYaxis()->SetTitleSize(0.06);
  hChi2X.GetXaxis()->SetLabelSize(0.06);
  hChi2X.GetYaxis()->SetLabelSize(0.06);

  hChi2Y.GetXaxis()->SetTitleSize(0.06);
  hChi2Y.GetYaxis()->SetTitleSize(0.06);
  hChi2Y.GetXaxis()->SetLabelSize(0.06);
  hChi2Y.GetYaxis()->SetLabelSize(0.06);


  TH1F hAngleBeforeChi2X = TH1F( Form("SinglePlaneAngleBeforeChi2CutX_ROC%i",plane_under_test), "SinglePlaneAngleBeforeChi2CutX", 100, -0.04, 0.04 );
  TH1F hAngleBeforeChi2Y = TH1F( Form("SinglePlaneAngleBeforeChi2CutY_ROC%i",plane_under_test), "SinglePlaneAngleBeforeChi2CutY", 100, -0.04, 0.04 );

  TH1F hAngleAfterChi2X = TH1F( Form("SinglePlaneAngleAfterChi2CutX_ROC%i",plane_under_test), "SinglePlaneAngleAfterChi2CutX", 100, -0.04, 0.04 );
  TH1F hAngleAfterChi2Y = TH1F( Form("SinglePlaneAngleAfterChi2CutY_ROC%i",plane_under_test), "SinglePlaneAngleAfterChi2CutY", 100, -0.04, 0.04 );


  double tz = FR->GetAlignment()->GetTZ(1, plane_under_test);

  // Event Loop
  for (int ievent = 0; FR->GetNextEvent() >= 0; ++ievent) {

    // print progress
    if (ievent % 10000 == 0) {
      std::cout << "Processing event: " << ievent << std::endl;
    }

    int i_slice = ievent/slice_size;
    if (i_slice==n_slices)
        i_slice--;

    // require exactly one track
    if (FR->NTracks() == 1){

      // Calculate the Angle of the tracks
      double slopeX = FR->Track(0)->fTVX / FR->Track(0)->fTVZ;
      double slopeY = FR->Track(0)->fTVY / FR->Track(0)->fTVZ;

      double angleX = atan(slopeX);
      double angleY = atan(slopeY);

      hAngleBeforeChi2X.Fill(angleX);
      hAngleBeforeChi2Y.Fill(angleY);

      hChi2.Fill( FR->Track(0)->Chi2());
      hChi2X.Fill( FR->Track(0)->Chi2X());
      hChi2Y.Fill( FR->Track(0)->Chi2Y());

      // Only accept reasonably central events
      if ((fabs(angleX) > 0.02) || (fabs(angleY) > 0.02))
        continue;


      // Look at the 90% quantile
      if (FR->Track(0)->Chi2X() > 6.25)
        continue;
      if (FR->Track(0)->Chi2Y() > 6.25)
        continue;

      hAngleAfterChi2X.Fill(angleX);
      hAngleAfterChi2Y.Fill(angleY);

      // Get the intersection of track and plane under test and fill
      // denominator histogram
      double tx = FR->Track(0)->TX( tz );
      double ty = FR->Track(0)->TY( tz );

      double lx = FR->GetAlignment()->TtoLX( tx, ty, 1, plane_under_test);
      double ly = FR->GetAlignment()->TtoLY( tx, ty, 1, plane_under_test);

      int px = FR->GetAlignment()->PXfromLX( lx );
      int py = FR->GetAlignment()->PYfromLY( ly );

      hOccupancyDenom.Fill( px, py );
      hOccupancyDenom_eventSlices[i_slice].Fill(px, py);

      PLTPlane* Plane = FR->Plane( plane_under_test );

      std::vector<float> delta_rs;

      for (uint16_t icl = 0; icl != Plane->NClusters(); icl++){

        float cl_px = Plane->Cluster(icl)->PX();
        float cl_py = Plane->Cluster(icl)->PY();

        float delta_px = px - cl_px;
        float delta_py = py - cl_py;

        delta_rs.push_back(sqrt(delta_px*delta_px + delta_py*delta_py));
      }

      if (DEBUG)
        std::cout << "TestPlaneEfficiency. Before FindILowestIndexAndValue." << std::endl;

      int closest_cluster_index = FindILowestIndexAndValue(delta_rs).first;

      if (DEBUG)
        std::cout << "TestPlaneEfficiency. After FindILowestIndexAndValue. closest_cluster_index = " << closest_cluster_index << std::endl;

      if (delta_rs.size() == 1)
        hDrSecondCluster.Fill(-1.);
      else if (delta_rs.size() >= 2)
        hDrSecondCluster.Fill(FindILowestIndexAndValue(delta_rs, 1).second);

      // Now look for a close hit in the plane under test

      int matched = 0;

      std::vector<float> charges_in_ell_1;
      std::vector<float> charges_in_ell_2;
      std::vector<float> charges_in_ell_3;
      std::vector<float> charges_in_ell_4;

      std::vector<int> adcs_in_ell_1;
      std::vector<int> adcs_in_ell_2;
      std::vector<int> adcs_in_ell_3;
      std::vector<int> adcs_in_ell_4;

      // Make sure there is at least one cluster
      if (closest_cluster_index != -1){

          // Determine here if the closest cluster is actually close enouigh
          // and fill cluster size
          float cluster_dtx = (tx - Plane->Cluster(closest_cluster_index)->TX());
          float cluster_dty = (ty - Plane->Cluster(closest_cluster_index)->TY());
          if (CheckEllipse(cluster_dtx, cluster_dty, max_dr_x, max_dr_y)){
            hClusterSize.Fill(px, py, Plane->Cluster(closest_cluster_index)->NHits());
          }

          if (DEBUG)
            std::cout << "TestPlaneEfficiency. Before Loop over hits" << std::endl;

          // loop over all hits in the cluster and check distance to intersection
          for (uint16_t ih = 0; ih != Plane->Cluster(closest_cluster_index)->NHits(); ih++){

                 if (DEBUG)
                   std::cout << "TestPlaneEfficiency. ih = " << ih << std::endl;

                 float dtx = (tx - Plane->Cluster(closest_cluster_index)->Hit(ih)->TX());
                 float dty = (ty - Plane->Cluster(closest_cluster_index)->Hit(ih)->TY());
                 float dtr = sqrt( dtx*dtx + dty*dty );

                 hdtx.Fill( dtx );
                 hdty.Fill( dty );
                 hdtr.Fill( dtr );

                 int adc = Plane->Cluster(closest_cluster_index)->Hit(ih)->ADC();
                 float charge = Plane->Cluster(closest_cluster_index)->Hit(ih)->Charge();

                 if (CheckEllipse(dtx, dty, max_dr_x, max_dr_y))
                    matched++;

                 // 1 Pixel Ellipse
                 if (CheckEllipse(dtx, dty, 0.015, 0.01)){
                   adcs_in_ell_1.push_back(adc);
                   charges_in_ell_1.push_back(charge);
                 }

                 // 2 Pixel Ellipse
                 if (CheckEllipse(dtx, dty, 0.03, 0.02)){
                   adcs_in_ell_2.push_back(adc);
                   charges_in_ell_2.push_back(charge);
                 }

                 // 3 Pixel Ellipse
                 if (CheckEllipse(dtx, dty, 0.045, 0.03)){
                   adcs_in_ell_3.push_back(adc);
                   charges_in_ell_3.push_back(charge);
                 }

                 // 4 Pixel Ellipse
                 if (CheckEllipse(dtx, dty, 0.06, 0.04)){
                   adcs_in_ell_4.push_back(adc);
                   charges_in_ell_4.push_back(charge);
                 }

          } // end of loop over hits

          hFractionContainted1.Fill(1. * charges_in_ell_1.size() / Plane->Cluster(closest_cluster_index)->NHits());
          hFractionContainted2.Fill(1. * charges_in_ell_2.size() / Plane->Cluster(closest_cluster_index)->NHits());
          hFractionContainted3.Fill(1. * charges_in_ell_3.size() / Plane->Cluster(closest_cluster_index)->NHits());
          hFractionContainted4.Fill(1. * charges_in_ell_4.size() / Plane->Cluster(closest_cluster_index)->NHits());

      } // End of having at least one valid cluster

      if (DEBUG)
        std::cout << "TestPlaneEfficiency. After Loop over hits" << std::endl;

      // if there was at least one match: fill denominator
      if (matched > 0){
         hOccupancyNum.Fill( px, py );
         hOccupancyNum_eventSlices[i_slice].Fill(px, py, 1);
      }

      // Sort Charge Vectors
      std::sort(charges_in_ell_1.begin(), charges_in_ell_1.end());
      std::sort(charges_in_ell_2.begin(), charges_in_ell_2.end());
      std::sort(charges_in_ell_3.begin(), charges_in_ell_3.end());
      std::sort(charges_in_ell_4.begin(), charges_in_ell_4.end());

      // Sort ADC Vectors
      std::sort(adcs_in_ell_1.begin(), adcs_in_ell_1.end());
      std::sort(adcs_in_ell_2.begin(), adcs_in_ell_2.end());
      std::sort(adcs_in_ell_3.begin(), adcs_in_ell_3.end());
      std::sort(adcs_in_ell_4.begin(), adcs_in_ell_4.end());

      // Fill Sum of Charges
      hSumCharge1.Fill(px, py, std::accumulate(charges_in_ell_1.begin(), charges_in_ell_1.end(), 0));
      hSumCharge2.Fill(px, py, std::accumulate(charges_in_ell_2.begin(), charges_in_ell_2.end(), 0));
      hSumCharge3.Fill(px, py, std::accumulate(charges_in_ell_3.begin(), charges_in_ell_3.end(), 0));
      hSumCharge4.Fill(px, py, std::accumulate(charges_in_ell_4.begin(), charges_in_ell_4.end(), 0));

      // Fill Highest Charge
      FillIth(&h1stCharge1, px, py, charges_in_ell_1, -1);
      FillIth(&h1stCharge2, px, py, charges_in_ell_2, -1);
      FillIth(&h1stCharge3, px, py, charges_in_ell_3, -1);
      FillIth(&h1stCharge4, px, py, charges_in_ell_4, -1);

      // Fill Highest ADC
      FillIth(&h1stCharge1ADC, px, py, adcs_in_ell_1, -1);
      FillIth(&h1stCharge2ADC, px, py, adcs_in_ell_2, -1);
      FillIth(&h1stCharge3ADC, px, py, adcs_in_ell_3, -1);
      FillIth(&h1stCharge4ADC, px, py, adcs_in_ell_4, -1);

      // Fill Second Highest Charge
      // do NOT fill if not available!
      FillIth(&h2ndCharge1, px, py, charges_in_ell_1, -2, false);
      FillIth(&h2ndCharge2, px, py, charges_in_ell_2, -2, false);
      FillIth(&h2ndCharge3, px, py, charges_in_ell_3, -2, false);
      FillIth(&h2ndCharge4, px, py, charges_in_ell_4, -2, false);

      // Fill Second Highest ADC
      // do NOT fill if not available!
      FillIth(&h2ndCharge1ADC, px, py, adcs_in_ell_1, -2, false);
      FillIth(&h2ndCharge2ADC, px, py, adcs_in_ell_2, -2, false);
      FillIth(&h2ndCharge3ADC, px, py, adcs_in_ell_3, -2, false);
      FillIth(&h2ndCharge4ADC, px, py, adcs_in_ell_4, -2, false);


      if (charges_in_ell_4.size() == 1){
	hCS1_1stCharge.Fill(charges_in_ell_4[0]);
	hCS1_SumCharge.Fill(charges_in_ell_4[0]);
      }

      if (charges_in_ell_4.size() == 2){
	  hCS2_1stCharge.Fill(charges_in_ell_4[1]);
	  hCS2_2ndCharge.Fill(charges_in_ell_4[0]);
	  hCS2_SumCharge.Fill(charges_in_ell_4[0]+charges_in_ell_4[1]);

	  if (rand() % 2 == 1)
	    hCS2_2D.Fill(charges_in_ell_4[0], charges_in_ell_4[1]);
	  else
	    hCS2_2D.Fill(charges_in_ell_4[1], charges_in_ell_4[0]);

      }

      if (charges_in_ell_4.size() == 3){
    	  hCS3_1stCharge.Fill(charges_in_ell_4[2]);
	  hCS3_2ndCharge.Fill(charges_in_ell_4[1]);
	  hCS3_3rdCharge.Fill(charges_in_ell_4[0]);
	  hCS3_SumCharge.Fill(charges_in_ell_4[0]+charges_in_ell_4[1]+charges_in_ell_4[2]);
      }

      if (charges_in_ell_4.size() == 4){
	hCS4_1stCharge.Fill(charges_in_ell_4[3]);
	hCS4_2ndCharge.Fill(charges_in_ell_4[2]);
	hCS4_3rdCharge.Fill(charges_in_ell_4[1]);
	hCS4_4thCharge.Fill(charges_in_ell_4[0]);
	hCS4_SumCharge.Fill(charges_in_ell_4[0]+charges_in_ell_4[1]+charges_in_ell_4[2]+charges_in_ell_4[3]);
      }





    } // end of having one track
  } // End of Event Loop




  // Remove masked areas from Occupancy Histograms
  const std::set<int> * pixelMask = FR->GetPixelMask();

  std::cout << "Got PixelMask: "<<pixelMask->size() <<std::endl;

  // Loop over all masked pixels
  for (std::set<int>::const_iterator ipix = pixelMask->begin();
       ipix != pixelMask->end();
       ipix++){

         // Decode the integer
         int roc  = (*ipix % 100000) / 10000;
         int col  = (*ipix % 10000 ) / 100;
         int row  = (*ipix % 100);

         // Make sure this concerns the plane under test
         if (roc == plane_under_test){

             // Convert pixel row/column to local coordinates
             // deltaR(local) should be == deltaR(telescope) (within a plane)
             float masked_lx = FR->GetAlignment()->PXtoLX( col);
             float masked_ly = FR->GetAlignment()->PYtoLY( row);

             //std::cout << col << " " << row << " " << masked_lx << " " << masked_ly << std::endl;

             // Loop over the TH2
             for (int ibin_x = 1; ibin_x != hOccupancyNum.GetNbinsX()+2; ibin_x++){
               for (int ibin_y = 1; ibin_y != hOccupancyNum.GetNbinsY()+2; ibin_y++){

                 // Get the bin-centers
                 int px =  hOccupancyNum.GetXaxis()->GetBinCenter( ibin_x );
                 int py =  hOccupancyNum.GetYaxis()->GetBinCenter( ibin_y );

                 float lx = FR->GetAlignment()->PXtoLX( px);
                 float ly = FR->GetAlignment()->PYtoLY( py);

                 //std::cout << px << " " << py << " " << lx << " " << ly;

                 // And check if they are within matching-distance of a masked pixel
                 float dtx = lx-masked_lx;
                 float dty = ly-masked_ly;

                 if (CheckEllipse(dtx, dty, max_dr_x, max_dr_y)){
                   // If yes: set numerator and denominator to zero
                   hOccupancyNum.SetBinContent( ibin_x, ibin_y, 0);
                   hOccupancyDenom.SetBinContent( ibin_x, ibin_y, 0);

                   for (int ibin_z = 1; ibin_z != hSumCharge1.GetNbinsZ()+2; ibin_z++){

                    hSumCharge1.SetBinContent(ibin_x, ibin_y, ibin_z, 0);
                    hSumCharge2.SetBinContent(ibin_x, ibin_y, ibin_z, 0);
                    hSumCharge3.SetBinContent(ibin_x, ibin_y, ibin_z, 0);
                    hSumCharge4.SetBinContent(ibin_x, ibin_y, ibin_z, 0);

                    h1stCharge1.SetBinContent(ibin_x, ibin_y, ibin_z, 0);
                    h1stCharge2.SetBinContent(ibin_x, ibin_y, ibin_z, 0);
                    h1stCharge3.SetBinContent(ibin_x, ibin_y, ibin_z, 0);
                    h1stCharge4.SetBinContent(ibin_x, ibin_y, ibin_z, 0);

                    h1stCharge1ADC.SetBinContent(ibin_x, ibin_y, ibin_z, 0);
                    h1stCharge2ADC.SetBinContent(ibin_x, ibin_y, ibin_z, 0);
                    h1stCharge3ADC.SetBinContent(ibin_x, ibin_y, ibin_z, 0);
                    h1stCharge4ADC.SetBinContent(ibin_x, ibin_y, ibin_z, 0);

                    h2ndCharge1.SetBinContent(ibin_x, ibin_y, ibin_z, 0);
                    h2ndCharge2.SetBinContent(ibin_x, ibin_y, ibin_z, 0);
                    h2ndCharge3.SetBinContent(ibin_x, ibin_y, ibin_z, 0);
                    h2ndCharge4.SetBinContent(ibin_x, ibin_y, ibin_z, 0);

                    h2ndCharge1ADC.SetBinContent(ibin_x, ibin_y, ibin_z, 0);
                    h2ndCharge2ADC.SetBinContent(ibin_x, ibin_y, ibin_z, 0);
                    h2ndCharge3ADC.SetBinContent(ibin_x, ibin_y, ibin_z, 0);
                    h2ndCharge4ADC.SetBinContent(ibin_x, ibin_y, ibin_z, 0);

                    hClusterSize.SetBinContent( ibin_x, ibin_y, ibin_z, 0);
                    for (int i=0; i != n_slices; i++){
                        hOccupancyNum_eventSlices[i].SetBinContent(ibin_x, ibin_y, ibin_z, 0);
                        hOccupancyDenom_eventSlices[i].SetBinContent(ibin_x, ibin_y, ibin_z, 0);
                    }

                   }

                 }

               }
             }
        }
   } // end loop over pixels


  // Prepare drawing
  TCanvas Can;
  Can.cd();
  out_f->cd();

  for (int i=0; i!= n_slices; i++){
    hOccupancyNum_eventSlices[i].Write();
    hOccupancyDenom_eventSlices[i].Write();
  }

  hOccupancyNum.SetMinimum(0);
  hOccupancyNum.SetAxisRange(12,38,"X");
  hOccupancyNum.SetAxisRange(39,80,"Y");
  hOccupancyNum.Draw("colz");
  hOccupancyNum.Write();


  hOccupancyDenom.Draw("colz");
  hOccupancyDenom.Write();
  Can.SaveAs( OutDir+TString(hOccupancyDenom.GetName()) + ".png");
  Can.SaveAs( OutDir+TString(hOccupancyDenom.GetName()) + ".pdf");

  hOccupancyDenom.SetMinimum(0);
  hOccupancyNum.SetAxisRange(12,38,"X");
  hOccupancyNum.SetAxisRange(39,80,"Y");
  hOccupancyDenom.SetAxisRange(12,38,"X");
  hOccupancyDenom.SetAxisRange(39,80,"Y");

  hOccupancyDenom.Draw("colz");
  hOccupancyDenom.Write();
  Can.SaveAs( OutDir+TString(hOccupancyDenom.GetName()) + ".png");
  Can.SaveAs( OutDir+TString(hOccupancyDenom.GetName()) + ".pdf");


  // Draw ratio of Occupancy histograms
  hOccupancyNum.Divide( &hOccupancyDenom );
  hOccupancyNum.SetMinimum(0);
  hOccupancyNum.SetMaximum(1.2);

  hOccupancyNum.Draw("colz");
  // Do not write the numerator-histo after division to the file
  Can.SaveAs( OutDir+TString(hOccupancyNum.GetName()) + ".png");
  Can.SaveAs( OutDir+TString(hOccupancyNum.GetName()) + ".pdf");

  hdtx.Draw();
  hdtx.Write();
  Can.SaveAs( OutDir+ TString(hdtx.GetName()) +".png");
  Can.SaveAs( OutDir+ TString(hdtx.GetName()) +".pdf");


  hdty.Draw();
  hdty.Write();
  Can.SaveAs(OutDir+ TString(hdty.GetName()) +".png");
  Can.SaveAs(OutDir+ TString(hdty.GetName()) +".pdf");

  hdtr.Draw();
  hdtr.Write();
  Can.SaveAs( OutDir+ TString(hdtr.GetName()) +".png");
  Can.SaveAs( OutDir+ TString(hdtr.GetName()) +".pdf");


  TF1 fun_chi2_6dof("chi2_6dof", "exp(-x/2.)*x*x/(4*16)");
  fun_chi2_6dof.SetRange(0.,50.);
  fun_chi2_6dof.SetNpx(1000);
  fun_chi2_6dof.Draw("SAME");

  hChi2.Scale(1/hChi2.Integral());
  hChi2.Draw();
  fun_chi2_6dof.Draw("SAME");
  hChi2.Write();
  Can.SaveAs( OutDir+ TString(hChi2.GetName()) +".png");
  Can.SaveAs( OutDir+ TString(hChi2.GetName()) +".pdf");

  hChi2X.SetTitle("");
  hChi2X.GetXaxis()->SetTitle("#chi^{2} (x)");
  hChi2X.GetYaxis()->SetTitle("A.U.");
  hChi2X.Scale( 1/ hChi2X.Integral());
  hChi2X.SetAxisRange(0, 0.07,"Y");
  hChi2X.SetLineWidth(2);
  hChi2X.Draw("hist");


  TF1 fun_chi2_3dof("chi2_3dof", "exp(-x/2.)*sqrt(x)/(5*sqrt(2*3.1415))");
  fun_chi2_3dof.SetRange(0.,20.);
  fun_chi2_3dof.SetNpx(1000);
  fun_chi2_3dof.SetLineWidth(2);
  fun_chi2_3dof.Draw("SAME");


  hChi2X.Write();
  Can.SaveAs( OutDir+ TString(hChi2X.GetName()) +".png");
  Can.SaveAs( OutDir+ TString(hChi2X.GetName()) +".pdf");


  hChi2Y.SetTitle("");
  hChi2Y.GetXaxis()->SetTitle("#chi^{2} (y)");
  hChi2Y.GetYaxis()->SetTitle("A.U.");
  hChi2Y.Scale(1/hChi2Y.Integral());
  hChi2Y.SetAxisRange(0, 0.08,"Y");
  hChi2Y.SetLineWidth(2);
  hChi2Y.Draw();
  fun_chi2_3dof.Draw("SAME");
  hChi2Y.Write();
  Can.SaveAs( OutDir+ TString(hChi2Y.GetName()) +".png");
  Can.SaveAs( OutDir+ TString(hChi2Y.GetName()) +".pdf");



  hCS1_1stCharge.Write();
  hCS2_1stCharge.Write();
  hCS2_2ndCharge.Write();
  hCS3_1stCharge.Write();
  hCS3_2ndCharge.Write();
  hCS3_3rdCharge.Write();
  hCS4_1stCharge.Write();
  hCS4_2ndCharge.Write();
  hCS4_3rdCharge.Write();
  hCS4_4thCharge.Write();


  hCS1_SumCharge.Write();
  hCS2_SumCharge.Write();
  hCS3_SumCharge.Write();
  hCS4_SumCharge.Write();

  hCS2_2D.Write();

  std::vector <TH3*> hs_mean_sum_charge;
  hs_mean_sum_charge.push_back( &hSumCharge1 );
  hs_mean_sum_charge.push_back( &hSumCharge2 );
  hs_mean_sum_charge.push_back( &hSumCharge3 );
  hs_mean_sum_charge.push_back( &hSumCharge4 );
  Write1DCharge(hs_mean_sum_charge, &Can, OutDir);

  std::vector <TH3*> hs_mean_1st_charge;
  hs_mean_1st_charge.push_back( &h1stCharge1 );
  hs_mean_1st_charge.push_back( &h1stCharge2 );
  hs_mean_1st_charge.push_back( &h1stCharge3 );
  hs_mean_1st_charge.push_back( &h1stCharge4 );
  Write1DCharge(hs_mean_1st_charge, &Can, OutDir);

  std::vector <TH3*> hs_mean_2nd_charge;
  hs_mean_2nd_charge.push_back( &h2ndCharge1 );
  hs_mean_2nd_charge.push_back( &h2ndCharge2 );
  hs_mean_2nd_charge.push_back( &h2ndCharge3 );
  hs_mean_2nd_charge.push_back( &h2ndCharge4 );
  Write1DCharge(hs_mean_2nd_charge, &Can, OutDir);

  std::vector <TH3*> hs_mean_1st_charge_adc;
  hs_mean_1st_charge_adc.push_back( &h1stCharge1ADC );
  hs_mean_1st_charge_adc.push_back( &h1stCharge2ADC );
  hs_mean_1st_charge_adc.push_back( &h1stCharge3ADC );
  hs_mean_1st_charge_adc.push_back( &h1stCharge4ADC );
  Write1DCharge(hs_mean_1st_charge_adc, &Can, OutDir);

  std::vector <TH3*> hs_mean_2nd_charge_adc;
  hs_mean_2nd_charge_adc.push_back( &h2ndCharge1ADC );
  hs_mean_2nd_charge_adc.push_back( &h2ndCharge2ADC );
  hs_mean_2nd_charge_adc.push_back( &h2ndCharge3ADC );
  hs_mean_2nd_charge_adc.push_back( &h2ndCharge4ADC );
  Write1DCharge(hs_mean_2nd_charge_adc, &Can, OutDir);


  float maxz;
  if (plane_under_test==1)
    maxz = 30000;
  if (plane_under_test==2)
    maxz = 30000;
  if (plane_under_test==3)
    maxz = 30000;
  if (plane_under_test==4)
    maxz = 50000;

  Write2DCharge( &hSumCharge1, &Can, maxz, OutDir);
  Write2DCharge( &hSumCharge2, &Can, maxz, OutDir);
  Write2DCharge( &hSumCharge3, &Can, maxz, OutDir);
  Write2DCharge( &hSumCharge4, &Can, maxz, OutDir);

  hSumCharge2.Write();
  hSumCharge4.Write();

  Write2DCharge( &h1stCharge1, &Can, maxz, OutDir);
  Write2DCharge( &h1stCharge2, &Can, maxz, OutDir);
  Write2DCharge( &h1stCharge3, &Can, maxz, OutDir);
  Write2DCharge( &h1stCharge4, &Can, maxz, OutDir);

  h1stCharge2.Write();
  h1stCharge4.Write();

  Write2DCharge( &h2ndCharge1, &Can, maxz, OutDir);
  Write2DCharge( &h2ndCharge2, &Can, maxz, OutDir);
  Write2DCharge( &h2ndCharge3, &Can, maxz, OutDir);
  Write2DCharge( &h2ndCharge4, &Can, maxz, OutDir);

  h2ndCharge2.Write();
  h2ndCharge4.Write();

  Write2DCharge( &hClusterSize, &Can, 7, OutDir);
  hClusterSize.Write();

  Can.SetLogy(1);
  hDrSecondCluster.Draw();
  Can.SaveAs( OutDir+ TString(hDrSecondCluster.GetName()) +".png");
  Can.SetLogy(0);

  std::vector<TH1*> hs_fraction_contained;
  hs_fraction_contained.push_back(&hFractionContainted1);
  hs_fraction_contained.push_back(&hFractionContainted2);
  hs_fraction_contained.push_back(&hFractionContainted3);
  hs_fraction_contained.push_back(&hFractionContainted4);
  Write1DFraction(hs_fraction_contained, &Can, OutDir);

  Can.SaveAs( OutDir+ TString(hFractionContainted1.GetName()) +".png");

  WriteAngleHistograms(&hAngleBeforeChi2X,
                       &hAngleBeforeChi2Y,
                       &hAngleAfterChi2X,
                       &hAngleAfterChi2Y,
                       &Can,
                       OutDir);

  delete FR;

}


void WriteHTML(TString const & OutDir, int telescopeID)
{
  // This function to write the HTML output for a run

  // Make output dir
  if (gSystem->mkdir(OutDir, true) != 0) {
    std::cerr << "WARNING: either OutDir exists or it is un-mkdir-able: " << OutDir << std::endl;
  }

  TString FileName;
  if (OutDir.Length() == 0) {
    FileName = OutDir + "/index.html";
  } else {
    FileName = OutDir + "/index.html";
  }
  std::ofstream f(FileName.Data());
  if (!f.is_open()) {
    std::cerr << "ERROR: Cannot open HTML file: " << FileName << std::endl;
    return;
  }



  f << "<html><body>\n";

  // RUN SUMMARY
  f << "<h1>Run Summary: </h1>\n";

  int nplanes = GetNPlanes();

  // LEVELS
  //f << "<hr />\n";
  //f << "<h2>Levels</h2>" << std::endl;
  //for (int i = 0; i != nplanes; ++i) {
  //  f << Form("<a href=\"Levels_ROC%i.png\"><img width=\"150\" src=\"Levels_ROC%i.png\"></a>\n", i, i);
  //}
  //f << "<br>" << std::endl;

  // OCCUPANCY
  f << "<hr />\n";
  f << "<h2>Occupancy</h2>" << std::endl;
  f << "<a href=\"Occupancy_Coincidence.png\"><img width=\"900\" src=\"Occupancy_Coincidence.png\"></a>\n<br>" << std::endl;
  for (int i = 0; i != nplanes; ++i) {
    f << Form("<a href=\"Occupancy_ROC%i.png\"><img width=\"150\" src=\"Occupancy_ROC%i.png\"></a>\n", i, i);
  }
  f << "<br>" << std::endl;
  for (int i = 0; i != nplanes; ++i) {
    f << Form("<a href=\"Occupancy_ROC%i_1DZ.png\"><img width=\"150\" src=\"Occupancy_ROC%i_1DZ.png\"></a>\n", i, i);
  }
  f << "<br>" << std::endl;
  for (int i = 0; i != nplanes; ++i) {
    f << Form("<a href=\"Occupancy_ROC%i_Quantile.png\"><img width=\"150\" src=\"Occupancy_ROC%i_Quantile.png\"></a>\n", i, i);
  }
  f << "<br>" << std::endl;
  for (int i = 0; i != nplanes; ++i) {
    f << Form("<a href=\"Occupancy1DZ_ROC%i_Quantile.png\"><img width=\"150\" src=\"Occupancy1DZ_ROC%i_Quantile.png\"></a>\n", i, i);
  }
  f << "<br>" << std::endl;

  for (int i = 0; i != nplanes; ++i) {
    f << Form("<a href=\"Occupancy_ROC%i_3x3Efficiency.png\"><img width=\"150\" src=\"Occupancy_ROC%i_3x3Efficiency.png\"></a>\n", i, i);
  }
  f << "<br>" << std::endl;
  for (int i = 0; i != nplanes; ++i) {
    f << Form("<a href=\"Occupancy_ROC%i_3x3Efficiency_1DZ.png\"><img width=\"150\" src=\"Occupancy_ROC%i_3x3Efficiency_1DZ.png\"></a>\n", i, i);
  }

  f << "<br>" << std::endl;
  for (int i = 0; i != nplanes; ++i) {
    f << Form("<a href=\"NClusters_ROC%i.png\"><img width=\"150\" src=\"NClusters_ROC%i.png\"></a>\n", i, i);
  }
  f << "<br>" << std::endl;
  f << "<br>" << std::endl;
  for (int i = 0; i != nplanes; ++i) {
    f << Form("<a href=\"NHitsPerCluster_ROC%i.png\"><img width=\"150\" src=\"NHitsPerCluster_ROC%i.png\"></a>\n", i, i);
  }
  f << "<br>" << std::endl;

  // PULSE HEIGHT
  f << "<hr />\n";
  f << "<h2>Pulse Height</h2>" << std::endl;
  for (int i = 0; i != nplanes; ++i) {
    f << Form("<a href=\"PulseHeight_ROC%i.png\"><img width=\"150\" src=\"PulseHeight_ROC%i.png\"></a>\n", i, i);
  }
  f << "<br>\n";
  for (int i = 0; i != nplanes; ++i) {
    f << Form("<a href=\"PulseHeightLong_ROC%i.png\"><img width=\"150\" src=\"PulseHeightLong_ROC%i.png\"></a>\n", i, i);
  }
  f << "<br>\n";
  for (int i = 0; i != nplanes; ++i) {
    f << Form("<a href=\"PulseHeightTime_ROC%i.png\"><img width=\"150\" src=\"PulseHeightTime_ROC%i.png\"></a>\n", i, i);
  }
  f << "<br>\n";
  for (int i = 0; i != nplanes; ++i) {
    f << Form("<a href=\"PulseHeightAvg2D_ROC%i.png\"><img width=\"150\" src=\"PulseHeightAvg2D_ROC%i.png\"></a>\n", i, i);
  }
  f << "<br>\n";
  for (int i = 0; i != nplanes; ++i) {
    f << Form("<a href=\"OccupancyLowPH_ROC%i.png\"><img width=\"150\" src=\"OccupancyLowPH_ROC%i.png\"></a>\n", i, i);
  }
  f << "<br>\n";
  for (int i = 0; i != nplanes; ++i) {
    f << Form("<a href=\"OccupancyHighPH_ROC%i.png\"><img width=\"150\" src=\"OccupancyHighPH_ROC%i.png\"></a>\n", i, i);
  }

  // TRACKING
  f << "<h2>Tracking</h2>\n";
  f << "<a href=\"TrackSlopeX.png\"><img width=\"150\" src=\"TrackSlopeX.png\"></a>\n";
  f << "<a href=\"TrackSlopeY.png\"><img width=\"150\" src=\"TrackSlopeY.png\"></a>\n";
  f << "<br>\n";
  f << "<a href=\"Chi2.png\"><img width=\"150\" src=\"Chi2.png\"></a>\n";
  f << "<a href=\"Chi2X.png\"><img width=\"150\" src=\"Chi2X.png\"></a>\n";
  f << "<a href=\"Chi2Y.png\"><img width=\"150\" src=\"Chi2Y.png\"></a>\n";

  // OFFLINE
//  f << "<h2>Straight Tracks</h2>\n";
//
//  f << "<br>\n";
//  for (int i = 0; i != nplanes; ++i) {
//    f << Form("<a href=\"PulseHeightOffline_ROC%i.png\"><img width=\"150\" src=\"PulseHeightOffline_ROC%i.png\"></a>\n", i, i);
//  }
//  f << "<br>\n";
//

  // TRACK RESIDUALS
  f << "<h2>Track Residuals</h2>\n";

  f << "<br>" << std::endl;
  for (int i = 0; i != nplanes; i++)
    f << Form("<a href=\"Residual_ROC%i_X.png\"><img width=\"150\" src=\"Residual_ROC%i_X.png\"></a>\n", i, i);
  f << "<br>\n";

  for (int i = 0; i != nplanes; i++)
    f << Form("<a href=\"Residual_ROC%i_Y.png\"><img width=\"150\" src=\"Residual_ROC%i_Y.png\"></a>\n", i, i);
  f << "<br>\n";

    /** Signal Distribution */
    if (GetNSignals() > 0){

        uint8_t nSig = GetNSignals();
        f << "<hr />\n";
        f << "<h2>Signal Distribution</h2>\n";
        for (uint8_t iSig = 0; iSig != nSig; iSig++)
            f << Form("<a href=\"Signal_%i.png\"><img width=\"150\" src=\"Signal_%i.png\"></a>\n", iSig, iSig);
        f << "<br>\n";
    }

  // Single Plane Studies
  if ( (telescopeID==1) || (telescopeID==2)){
    f << "<h2>Single Plane Studies</h2>\n";
    f << "<br>" << std::endl;

    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"TracksPassing_ROC%i.png\"><img width=\"150\" src=\"TracksPassing_ROC%i.png\"></a>\n", i, i);
    f << "<br>\n";


    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"PlaneEfficiency_ROC%i.png\"><img width=\"150\" src=\"PlaneEfficiency_ROC%i.png\"></a>\n", i, i);
    f << "<br>\n";

    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"ClusterSize_ROC%i_profile.png\"><img width=\"150\" src=\"ClusterSize_ROC%i_profile.png\"></a>\n", i, i);
    f << "<br>\n";


    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"SumCharge_ROC%i.png\"><img width=\"150\" src=\"SumCharge_ROC%i.png\"></a>\n", i, i);
    f << "<br>\n";

    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"1stCharge_ROC%i.png\"><img width=\"150\" src=\"1stCharge_ROC%i.png\"></a>\n", i, i);
    f << "<br>\n";

    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"2ndCharge_ROC%i.png\"><img width=\"150\" src=\"2ndCharge_ROC%i.png\"></a>\n", i, i);
    f << "<br>\n";



    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"SumCharge_ROC%i_profile.png\"><img width=\"150\" src=\"SumCharge_ROC%i_profile.png\"></a>\n", i, i);
    f << "<br>\n";

    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"SumCharge2_ROC%i_profile.png\"><img width=\"150\" src=\"SumCharge2_ROC%i_profile.png\"></a>\n", i, i);
    f << "<br>\n";

    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"SumCharge3_ROC%i_profile.png\"><img width=\"150\" src=\"SumCharge3_ROC%i_profile.png\"></a>\n", i, i);
    f << "<br>\n";

    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"SumCharge3_ROC%i_profile.png\"><img width=\"150\" src=\"SumCharge4_ROC%i_profile.png\"></a>\n", i, i);
    f << "<br>\n";


    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"1stCharge_ROC%i_profile.png\"><img width=\"150\" src=\"1stCharge_ROC%i_profile.png\"></a>\n", i, i);
    f << "<br>\n";

    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"1stCharge2_ROC%i_profile.png\"><img width=\"150\" src=\"1stCharge2_ROC%i_profile.png\"></a>\n", i, i);
    f << "<br>\n";

    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"1stCharge3_ROC%i_profile.png\"><img width=\"150\" src=\"1stCharge3_ROC%i_profile.png\"></a>\n", i, i);
    f << "<br>\n";

    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"1stCharge3_ROC%i_profile.png\"><img width=\"150\" src=\"1stCharge4_ROC%i_profile.png\"></a>\n", i, i);
    f << "<br>\n";

    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"2ndCharge_ROC%i_profile.png\"><img width=\"150\" src=\"2ndCharge_ROC%i_profile.png\"></a>\n", i, i);
    f << "<br>\n";

    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"2ndCharge2_ROC%i_profile.png\"><img width=\"150\" src=\"2ndCharge2_ROC%i_profile.png\"></a>\n", i, i);
    f << "<br>\n";

    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"2ndCharge3_ROC%i_profile.png\"><img width=\"150\" src=\"2ndCharge3_ROC%i_profile.png\"></a>\n", i, i);
    f << "<br>\n";

    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"2ndCharge3_ROC%i_profile.png\"><img width=\"150\" src=\"2ndCharge4_ROC%i_profile.png\"></a>\n", i, i);
    f << "<br>\n";




    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"SinglePlaneTestChi2_ROC%i.png\"><img width=\"150\" src=\"SinglePlaneTestChi2_ROC%i.png\"></a>\n", i, i);
    f << "<br>\n";

    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"SinglePlaneTestChi2X_ROC%i.png\"><img width=\"150\" src=\"SinglePlaneTestChi2X_ROC%i.png\"></a>\n", i, i);
    f << "<br>\n";

    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"SinglePlaneTestChi2Y_ROC%i.png\"><img width=\"150\" src=\"SinglePlaneTestChi2Y_ROC%i.png\"></a>\n", i, i);
    f << "<br>\n";


    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"SinglePlaneTestDY_ROC%i.png\"><img width=\"150\" src=\"SinglePlaneTestDY_ROC%i.png\"></a>\n", i, i);
    f << "<br>\n";

    for (int i = 1; i != 5; i++)
      f << Form("<a href=\"SinglePlaneTestDR_ROC%i.png\"><img width=\"150\" src=\"SinglePlaneTestDR_ROC%i.png\"></a>\n", i, i);
    f << "<br>\n";
  }

  // EVENT DISPLAYS
  f << "<h2>Event Displays</h2>\n";

  f << "<br>" << std::endl;
  for (int irow = 0; irow != 4; irow++){
    for (int icol = 1; icol != nplanes; ++icol) {
      int i = irow*5+icol;
      f << Form("<a href=\"Tracks_Ev%i.png\"><img width=\"150\" src=\"Tracks_Ev%i.png\"></a>\n", i, i);
    }
    f << "<br>\n";
  }
  f << "<br>\n";


  f << "</body></html>";
  f.close();


}


void PrintUsage(const string & name) {
  cerr << "Usage: " << name << " <InFileName> <action> <telescopeID> ";
  cerr << "optional arguments: (<TrackMode>=0) (<EventsAlignment>=100000) (<IterAlignStep>=20) (<MaxAlignRes(cm)>=0.00001) (<MaxAlignAngle(rad)>=0.001) (<SilDUT>=-1)" << endl;
  cerr << "action:\n  0: analysis\n  1: alignment\n  2: residuals" << endl;
  cerr << "TrackMode:\n  0: AllPlanes\n  1: OnlyTelescope" << endl;
  cerr << "EventsAlignment:\n  0: Use ALL events in file\n  <n>: Use only the first \"n\" events in the provided file." << endl;
  cerr << "SilDUT:\n  -1: No Silicon DUT, only diamonds\n  <i>: Roc position \"i\" where the Silicon DUT is" << endl;
}


int main (int argc, char* argv[]) {

  const uint16_t max_args = 11;
  if (argc <= 3 or argc >= max_args) {
    tel::critical("Wrong arguments; Must supply at least 3 arguments and no more than 9: ");
    PrintUsage(argv[0]);
    return 1;
  }
  gInterpreter->GenerateDictionary("vector<vector<float> >;vector<vector<UShort_t> >", "vector"); // add root dicts for vector<vector> >
  gROOT->ProcessLine("#include <vector>");

  /** There three usage modes: analysis, alignment and residuals
      analysis: uses alignment and residuals for the given telescope to perform global and single plane studies
      alignment: starts with all alignment constants zero and does several iterations to minimize the residuals. All planes are shifted in x and y and rotated
        around the z-axis. Residual plots of the last iteration are saved.
      residuals: tries to find the correct residuals for tracking
      action:
        0: Analysis
        1: Alignment
        2: Residuals */
  auto action = stoi(argv[2]);
  auto telescope_id = stoi(argv[3]); /** see data/alignments.txt file */
  /** Tracking only on the telescope (only for digital telescope):
      0: Use All planes (default until September 2016.
      1: Use only the first 4 planes for tracking (telescope planes) */
  auto track_only_telescope = argc >= 5 and bool(stoi(argv[4]));

  if (action > 3) {
    tel::critical("Wrong action argument: " + to_string(action));
    return 2;
  }
  vector<string> action_str = {" (Analysis)", " (Alignment)", " (Residuals)"};
  cout << "Action = " << int(action) << action_str.at(action) << endl;
  cout << "TelescopeID = " << int(telescope_id) << endl;
  cout << "Track only analogue telescope: " << track_only_telescope << endl << endl;

  /** read config */
  if (tel::Config::Read(telescope_id) == 0) { return 3; }

  /** optional settings */
  AlignSettings AS = ReadAlignSettings(vector<string>(argv, argv + argc), action_str.size());

  const string in_file_name = argv[1];
  const string run_number = tel::trim(tel::trim(tel::split(in_file_name, '/').back(), "estro0"), ".");

  ValidateDirectories(run_number);

  /** Open a ROOT file to store histograms in. */
  TFile out_f(Form("%s/plots/%s/histos.root", GetDir().c_str(), run_number.c_str()), "recreate");

  if (action == 1) { /** ALIGNMENT */
    Alignment(in_file_name, run_number, telescope_id, track_only_telescope, AS.n_iterations_, AS.res_thresh_, AS.angle_thresh_, AS.max_events_, AS.sil_roc_);
  } else if (action==2) { /** RESIDUAL CALCULATION */
    FindPlaneErrors(in_file_name, run_number, telescope_id);
  } else { /** ANALYSIS */
    PLTAnalysis Analysis(in_file_name, &out_f, run_number, telescope_id, bool(track_only_telescope), AS.max_events_);
    Analysis.EventLoop();
    Analysis.FinishAnalysis();
  }

  return 0;
}
