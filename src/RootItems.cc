#include "RootItems.h"

using namespace std;

/** ============================
 CONSTRUCTOR
 =================================*/
RootItems::RootItems(uint8_t telescopeID, TString const RunNumber):
    nRoc(GetNumberOfROCS(telescopeID)),
    PlotsDir("plots/"),
    OutDir(PlotsDir + RunNumber + "/"),
    HistColors {1, 4, 28, 2 } {

    /** canvases */
    c1 = new TCanvas;
    c2 = new TCanvas("CoincidenceMap", "CoincidenceMap", 1200, 400);

    /** tracking */
    hTrackSlopeX = new TH1F("TrackSlopeX", "TrackSlopeX", 50, -0.05, 0.05);
    hTrackSlopeY = new TH1F("TrackSlopeY", "TrackSlopeY", 50, -0.05, 0.05);
    fGauss = new TF1("fGauss", "gaus", -0.05, 0.05);
    lFitGauss = new TLegend(0.7,0.65,0.88,0.85);

    /** occupancy */
    hOccupancy = FillVectorTH2F(hOccupancy, "Occupancy_ROC%i");
    hOccupancyLowPH = FillVectorTH2F(hOccupancyLowPH, "OccupancyLowPH_ROC%i");
    hOccupancyHighPH = FillVectorTH2F(hOccupancyHighPH, "OccupancyHighPH_ROC%i");

    /** cluster */
    hNHitsPerCluster = FillVectorTH1F(hNHitsPerCluster, "NHitsPerCluster_ROC%i");
    hNClusters = FillVectorTH1F(hNClusters, "NClusters_ROC%i");

    /** pulse height */
    hPulseHeight = FillVectorPH(hPulseHeight, "", 50000);
    hPulseHeightLong = FillVectorPH(hPulseHeightLong, "Long", 300000);
    hPulseHeightOffline = FillVectorPH(hPulseHeightOffline, "Offline", 50000);
    lPulseHeight = new TLegend(0.77, 0.7, 0.9, 0.88, "");
    lPHMean = new TLegend(0.77, 0.45, 0.9, 0.66, "Mean:");

    /** coincidence map */
    hCoincidenceMap = new TH1F("CoincidenceMap", "CoincidenceMap", pow(2, nRoc), 0, pow(2, nRoc));
}
RootItems::~RootItems() { }


/** ============================
 AUXILIARY FUNCTIONS
 =================================*/
void RootItems::FitSlope(TH1F * histo){

    fGauss->SetLineWidth(2);
    histo->Fit(fGauss, "Q");
}
void RootItems::LegendSlope(TH1F * histo) {

    lFitGauss->Clear();
    lFitGauss->SetTextSize(0.03);
    lFitGauss->AddEntry(histo,"Track Slope X","l");
    lFitGauss->AddEntry(fGauss,"Global Fit","l");
    TString s1 = "Mean: ";
    TString s2 = "Max: ";
    TString s3 = "Sigma: ";
    lFitGauss->AddEntry("Mean", s1 + TString::Format("%1.3f", fGauss->GetParameter(1)), "");
    lFitGauss->AddEntry("const", s2 + TString::Format("%5.1f", fGauss->GetParameter(0)), "");
    lFitGauss->AddEntry("const", s3 + TString::Format("%1.3f", fGauss->GetParameter(2)), "");
    lFitGauss->Draw("same");
}
vector<TH2F*> RootItems::FillVectorTH2F(vector<TH2F*> histo, const char * name) {
    for (uint16_t iroc = 0; iroc != nRoc; ++iroc){
        TH2F * hist = new TH2F(Form(name, iroc), Form(name, iroc), 52, 0, 52, 80, 0, 80);
        histo.push_back(hist);
    }
    return histo;
}
vector<TH1F*> RootItems::FillVectorTH1F(vector<TH1F*> histo, const char * name) {
    for (uint16_t iroc = 0; iroc != nRoc; ++iroc){
        TH1F * hist = new TH1F(Form(name, iroc), Form(name, iroc), 10, 0, 10);
        histo.push_back(hist);
    }
    return histo;
}
std::vector<vector<TH1F*> > RootItems::FillVectorPH(std::vector<vector<TH1F*> > histVec, TString name, uint32_t maxPH){
    const uint8_t minPH(0), nBins(50);
    TString base = "PulseHeight" + name + "_ROC%i_";
    TString names[4] = {base + "All", base + "NPix1", base + "NPix2", base + "NPix3Plus"};
    for (uint8_t iroc(0); iroc != nRoc; ++iroc){
        histVec.resize(iroc + 1);
        for (uint8_t iMode(0); iMode != 4; iMode++){
            TString Name = TString::Format(names[iMode], iroc);
            TH1F * hist = new TH1F(Name, Name, nBins, minPH, maxPH);
            histVec[iroc].push_back(hist);
        }
    }
    return histVec;
}
void RootItems::FormatPHHisto(std::vector<vector<TH1F*> > histVec){

    for (uint8_t iroc(0); iroc != nRoc; ++iroc)
        for (uint8_t iMode(0); iMode != 4; iMode++){
            histVec[iroc][iMode]->SetXTitle("Charge (electrons)");
            histVec[iroc][iMode]->SetYTitle("Number of Clusters");
            histVec[iroc][iMode]->SetLineColor(HistColors[iMode]);
    }
}
void RootItems::FormatLegendPH(){

    lPulseHeight->SetFillColor(4000);
    lPulseHeight->SetFillStyle(4000);
    lPulseHeight->SetBorderSize(0);
    lPulseHeight->SetTextAlign(11);
    lPHMean->SetTextAlign(11);
    lPHMean->SetFillStyle(0);
    lPHMean->SetBorderSize(0);
}
void RootItems::FillLegendsPH(uint8_t iroc, std::vector<vector<TH1F*> > histVec){

    TString names1[4] = {"All", "1 Pix", "2 Pix", "3+ Pix"};
    TString names2[4] = {"PH0PMean", "PH1PMean", "PH2PMean", "PH3PMean"};
    for (uint8_t iMode(0); iMode != 4; iMode++){
        lPulseHeight->AddEntry(histVec[iroc][iMode], names1[iMode], "l");
        lPHMean->AddEntry(names2[iMode], TString::Format("%8.0f", histVec[iroc][iMode]->GetMean()), "")->SetTextColor(HistColors[iMode]);
    }
}
void RootItems::DrawSavePH(uint8_t iroc, std::vector<vector<TH1F*> > histVec, TString title, TString saveName){

    c1->cd();
    histVec[iroc][0]->SetTitle( TString::Format(title, iroc));
    histVec[iroc][0]->Draw("hist");
    for (uint8_t i(1); i != 4; i++) histVec[iroc][i]->Draw("samehist");
    lPulseHeight->Draw("same");
    lPHMean->Draw("same");
    c1->SaveAs(OutDir+TString::Format(saveName, iroc));
    for (uint8_t i(0); i != 4; i++) histVec[iroc][i]->Write();
}
void RootItems::ClearLegendsPH(){
    lPulseHeight->Clear();
    lPHMean->Clear();
    lPHMean->SetHeader("Mean:");
}
void RootItems::DrawSaveTH1F(std::vector<TH1F*> histo, uint8_t iroc, TCanvas & c, const char * xTit, const char * yTit){
    c.cd();
    histo[iroc]->SetMinimum(0);
    histo[iroc]->SetXTitle(xTit);
    histo[iroc]->SetYTitle(yTit);
    histo[iroc]->Draw("hist");
    c.SaveAs(OutDir+TString(histo[iroc]->GetName()) + ".gif");
}
void RootItems::PrepCoincidenceHisto(){
    hCoincidenceMap->SetFillColor(40);
    hCoincidenceMap->GetYaxis()->SetTitle("Number of Hits");
    hCoincidenceMap->GetYaxis()->SetTitleOffset(0.5);
    hCoincidenceMap->GetYaxis()->SetTitleSize(0.05);
    hCoincidenceMap->GetXaxis()->SetTitle("Coincidences");
    hCoincidenceMap->GetXaxis()->SetTitleSize(0.05);
    hCoincidenceMap->GetXaxis()->SetLabelOffset(99);
}
void RootItems::DrawSaveCoincidence(){
    c2->cd();
    c2->SetLogy(1);
    hCoincidenceMap->Draw("");
    vector<string> bin;
    for (uint16_t i(0); i < pow(2, nRoc); i++){
        bin.push_back("");
        int z = i;
        for (int16_t j(nRoc - 1); j >= 0; j--){
            int x = pow(2, j);
            int y = z / x;
            if (y == 1) z -= x;
            bin[i].push_back(y +'0');
        }
    }
    TText Labels;
    Labels.SetTextAngle(0);
    Labels.SetTextSize(0.04);
    Labels.SetTextAlign(22);
    Float_t x, y(0);
    for (uint8_t iBins(0); iBins < pow(2, nRoc); iBins++) {
        x = hCoincidenceMap->GetXaxis()->GetBinCenter(iBins + 1);
        Labels.DrawText(x, y, bin[iBins].c_str());}
    c2->SaveAs(OutDir+"Occupancy_Coincidence.gif");
}
