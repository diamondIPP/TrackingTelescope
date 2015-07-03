#include "RootItems.h"

using namespace std;

/** constructor */
RootItems::RootItems() {
    hTrackSlopeX = new TH1F("TrackSlopeX", "TrackSlopeX", 50, -0.05, 0.05);
    hTrackSlopeY = new TH1F("TrackSlopeY", "TrackSlopeY", 50, -0.05, 0.05);
    fGauss = new TF1("fGauss", "gaus", -0.05, 0.05);
    lFitGauss = new TLegend(0.7,0.65,0.88,0.85);
}
RootItems::~RootItems() { }

/** get-functions*/
TH1F * RootItems::TrackSlopeX(){ return hTrackSlopeX; }
TH1F * RootItems::TrackSlopeY(){ return hTrackSlopeY; }
TLegend * RootItems::FitGauss(){ return lFitGauss; }

/** helper functions */

void RootItems::FitSlope(){

    fGauss->SetLineWidth(2);
    hTrackSlopeX->Fit(fGauss, "U");
    hTrackSlopeY->Fit(fGauss, "U");
}

void RootItems::LegendSlope(TH1F * histo) {

    lFitGauss->Clear();
    lFitGauss->SetTextSize(0.03);
    lFitGauss->AddEntry(histo,"Track Slope X","l");
    lFitGauss->AddEntry(fGauss,"Global Fit","l");
    lFitGauss->Draw("same");
}

