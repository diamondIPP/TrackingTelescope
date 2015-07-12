#include "RootItems.h"

using namespace std;

/** ============================
 CONSTRUCTOR
 =================================*/
RootItems::RootItems(uint8_t telescopeID, TString const RunNumber):
    nRoc(GetNumberOfROCS(telescopeID)),
    PlotsDir("plots/"),
    OutDir(PlotsDir + RunNumber + "/"),
    FileType(".gif"),
    HistColors {1, 4, 28, 2 },
    maxChi2(20) {

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
    hOccupancy1DZ.resize(nRoc);
    h3x3.resize(nRoc);
    h3x31DZ.resize(nRoc);
    hOccupancyLowPH = FillVectorTH2F(hOccupancyLowPH, "OccupancyLowPH_ROC%i");
    hOccupancyHighPH = FillVectorTH2F(hOccupancyHighPH, "OccupancyHighPH_ROC%i");

    /** cluster */
    hNHitsPerCluster = FillVectorTH1F(hNHitsPerCluster, "NHitsPerCluster_ROC%i");
    hNClusters = FillVectorTH1F(hNClusters, "NClusters_ROC%i");

    /** pulse height */
    hPulseHeight = FillVectorPH(hPulseHeight, "", 50000);
    FormatPHHisto(hPulseHeight);
    hPulseHeightLong = FillVectorPH(hPulseHeightLong, "Long", 300000);
    FormatPHHisto(hPulseHeightLong);
    hPulseHeightOffline = FillVectorPH(hPulseHeightOffline, "Offline", 50000);
    FormatPHHisto(hPulseHeightOffline);
    lPulseHeight = new TLegend(0.77, 0.7, 0.9, 0.88, "");
    lPHMean = new TLegend(0.77, 0.45, 0.9, 0.66, "Mean:");
    lRatio = new TLegend(0.77, 0.25, 0.9, 0.45, "Ratio [%]:");
    FormatLegendsPH();
    AllocateArrAvPH();
    gAvgPH = FillVecAvPH(gAvgPH);
    hPulseHeightAvg2D = FillVectorTH2F(hPulseHeightAvg2D, "PulseHeightAvg2D_ROC%i");

    /** coincidence map */
    hCoincidenceMap = new TH1F("CoincidenceMap", "CoincidenceMap", pow(2, nRoc), 0, pow(2, nRoc));

    /** chi2 */
    hChi2 = new TH1F("Chi2", "Chi2", 240, 0., maxChi2);
    hChi2X = new TH1F("Chi2X", "Chi2X", 240, 0., maxChi2);
    hChi2Y = new TH1F("Chi2Y", "Chi2Y", 240, 0., maxChi2);

    /** residuals */
    hResidual = FillVecResidual(hResidual, "Residual_ROC%i",  100, -.15, .15, 100, -.15, .15);
    hResidualXdY = FillVecResidual(hResidualXdY, "ResidualXdY_ROC%i", 200, -1, 1, 100, -.5, .5);
    hResidualYdX = FillVecResidual(hResidualYdX, "ResidualYdX_ROC%i", 200, -1, 1, 100, -.5, .5);


}
RootItems::~RootItems() {

    delete c1; delete c2;
    delete hTrackSlopeX; delete hTrackSlopeY; delete fGauss; delete lFitGauss;
    for (uint8_t iRoc = 0; iRoc != nRoc; iRoc++){
        delete hOccupancy[iRoc]; delete hOccupancy1DZ[iRoc]; delete hOccupancyLowPH[iRoc]; delete hOccupancyHighPH[iRoc];
        delete h3x3[iRoc]; delete h3x31DZ[iRoc];
    }
}


/** ============================
 MAIN FUNCTIONS
 =================================*/
 void RootItems::SaveAllHistos(){

    for (int iRoc = 0; iRoc != nRoc; ++iRoc) {

    /** occupancy */
    DrawSaveOccupancy(iRoc, Occupancy() );
    Occupancy()[iRoc]->Write();
    setOccupancy1DZ(PLTU::HistFrom2D(Occupancy()[iRoc]), iRoc);
    DrawSaveOccupancy1DZ(iRoc);
    Occupancy1DZ()[iRoc]->Write();
    DrawSaveOccupancyQuantile(iRoc);
    setOccupancy1DZ(PLTU::HistFrom2D(Occupancy()[iRoc], 0, QuantileValue()[0],
        TString::Format("Occupancy1DZ_ROC%i_Quantile", iRoc), 20), iRoc);
    DrawSaveOccupancy1DZ(iRoc);
    /** 3x3 efficiency */
    set3x3(PLTU::Get3x3EfficiencyHist(*Occupancy()[iRoc], 0, 51, 0, 79), iRoc);
    DrawSave3x3(iRoc);
    set3x31DZ(PLTU::HistFrom2D(Eff3x3()[iRoc], "", 50), iRoc);
    DrawSave3x31DZ(iRoc);
    /** low + high PH */
    DrawSaveOccupancy(iRoc, OccupancyLowPH() );
    DrawSaveOccupancy(iRoc, OccupancyHighPH() );

    /** clusters per event */
    DrawSaveTH1F(nClusters(), iRoc, "Number of clusters per event", "Events");

    /** hits per cluster */
    DrawSaveTH1F(nHitsPerCluster(), iRoc, "Number of hits per cluster", "Number of Clusters");

    /** pulse heights */
    /** standard */
    ClearLegendsPH();
    FillLegendsPH(iRoc, PulseHeight());
    DrawSavePH(iRoc, PulseHeight(), "Pulse Height ROC%i", "PulseHeight_ROC%i.gif");
    /** offline */
    ClearLegendsPH();
    FillLegendsPH(iRoc, PulseHeightOffline());
    DrawSavePH(iRoc, PulseHeightOffline(), "Pulse Height Offline ROC%i", "PulseHeightOffline_ROC%i.gif");
    /** long */
    ClearLegendsPH();
    FillLegendsPH(iRoc, PulseHeightLong());
    DrawSavePH(iRoc, PulseHeightLong(), "Pulse Height Long ROC%i", "PulseHeightLong_ROC%i.gif");
    /** average pulse height */
    DrawSaveAvPH(iRoc);
    FillAvPH2D(iRoc);
    DrawSaveAvPH2D(iRoc);

    /** residuals */
    DrawSaveResidual(iRoc, Residual());
    DrawSaveResidual(iRoc, ResidualXdY());
    DrawSaveResidual(iRoc, ResidualYdX());
    DrawSaveResidualProj(iRoc, Residual(), "X");
    DrawSaveResidualProj(iRoc, Residual(), "Y");

  } // end of loop over ROCs

    /** draw and save coincidence map */
    PrepCoincidenceHisto();
    DrawSaveCoincidence();

    /** draw tracking slopes and chi2 */
    DrawSaveTrackSlope(TrackSlopeX() );
    TrackSlopeX()->Write();
    DrawSaveTrackSlope(TrackSlopeY() );
    TrackSlopeY()->Write();
    DrawSaveChi2(Chi2(), "Chi2");
    DrawSaveChi2(Chi2X(), "Chi2X");
    DrawSaveChi2(Chi2Y(), "Chi2Y");
 }


/** ============================
 AUXILIARY FUNCTIONS
 =================================*/
void RootItems::FitSlope(TH1F * histo){

    fGauss->SetLineWidth(2);
    histo->Fit(fGauss, "Q");
//    histo->SetStats(true);
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
        TH2F * hist = new TH2F(Form(name, iroc), Form(name, iroc), PLTU::NCOL, PLTU::FIRSTCOL, PLTU::LASTCOL, PLTU::NROW, PLTU::FIRSTROW, PLTU::LASTROW);
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
vector<vector<TH1F*> > RootItems::FillVectorPH(vector<vector<TH1F*> > histVec, TString name, uint32_t maxPH){
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
void RootItems::FillAvPH2D(uint8_t iroc){

    for (uint8_t iCol = 0; iCol != PLTU::NCOL; ++iCol)
        for (uint8_t iRow = 0; iRow != PLTU::NROW; ++iRow)
            hPulseHeightAvg2D[iroc]->SetBinContent(iCol+1, iRow+1, dAvgPH2D[iroc][iCol][iRow]);
}
void RootItems::FormatPHHisto(std::vector<vector<TH1F*> > histVec){

    for (uint8_t iroc(0); iroc != nRoc; ++iroc)
        for (uint8_t iMode(0); iMode != 4; iMode++){
            histVec[iroc][iMode]->SetXTitle("Charge (electrons)");
            histVec[iroc][iMode]->SetYTitle("Number of Clusters");
            histVec[iroc][iMode]->SetLineColor(HistColors[iMode]);
    }
}
void RootItems::FormatLegendsPH(){

    lPulseHeight->SetFillColor(4000);
    lPulseHeight->SetFillStyle(4000);
    lPulseHeight->SetBorderSize(0);
    lPulseHeight->SetTextAlign(11);
    lPHMean->SetTextAlign(11);
    lPHMean->SetFillStyle(4000);
    lPHMean->SetBorderSize(0);
    lRatio->SetTextAlign(11);
    lRatio->SetFillStyle(4000);
    lRatio->SetBorderSize(0);
}
void RootItems::FillLegendsPH(uint8_t iroc, std::vector<vector<TH1F*> > histVec){

    TString names1[4] = {"All", "1 Pix", "2 Pix", "3+ Pix"};
    TString names2[4] = {"PH0PMean", "PH1PMean", "PH2PMean", "PH3PMean"};
    for (uint8_t iMode(0); iMode != 4; iMode++){
        lPulseHeight->AddEntry(histVec[iroc][iMode], names1[iMode], "l");
        lPHMean->AddEntry(names2[iMode], TString::Format("%8.0f", histVec[iroc][iMode]->GetMean()), "")->SetTextColor(HistColors[iMode]);
    }
    Float_t a = hPulseHeight[iroc][1]->GetEntries();
    Float_t b = hPulseHeight[iroc][2]->GetEntries();
    Float_t c = hPulseHeight[iroc][3]->GetEntries();
    lRatio->AddEntry("2/1", "2/1: " + TString::Format("%02.2f", b / a * 100), "");
    lRatio->AddEntry("3/1", "3/1: " + TString::Format("%02.2f", c / a * 100), "");
    lRatio->AddEntry("3/2", "3/2: " + TString::Format("%02.2f", c / b * 100), "");
}
void RootItems::DrawSavePH(uint8_t iroc, std::vector<vector<TH1F*> > histVec, TString title, TString saveName){

    c1->cd();
    histVec[iroc][0]->SetTitle( TString::Format(title, iroc));
    histVec[iroc][0]->Draw("hist");
    for (uint8_t i(1); i != 4; i++) histVec[iroc][i]->Draw("samehist");
    lPulseHeight->Draw("same");
    lPHMean->Draw("same");
    lRatio->Draw("same");
    c1->SaveAs(OutDir+TString::Format(saveName, iroc));
    for (uint8_t i(0); i != 4; i++) histVec[iroc][i]->Write();
}
void RootItems::ClearLegendsPH(){
    lPulseHeight->Clear();
    lPHMean->Clear();
    lRatio->Clear();
    lPHMean->SetHeader("Mean:");
    lRatio->SetHeader("Ratio: [%]");
}
void RootItems::DrawSaveTH1F(std::vector<TH1F*> histo, uint8_t iroc, const char * xTit, const char * yTit){

    c1->cd();
    histo[iroc]->SetMinimum(0);
    histo[iroc]->SetXTitle(xTit);
    histo[iroc]->SetYTitle(yTit);
    histo[iroc]->Draw("hist");
    c1->SaveAs(OutDir+TString(histo[iroc]->GetName()) + ".gif");
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
            bin[i].push_back(y + '0');
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
    c2->SaveAs(OutDir + "Occupancy_Coincidence.gif");
}
void RootItems::DrawSaveChi2(TH1F * histo, TString saveName){
    c1->cd();
//    hChi2X.Scale( 1/hChi2X.Integral());
    gStyle->SetOptStat(0);
    histo->Draw("hist");
    c1->SaveAs(OutDir + saveName + ".gif");
    histo->Write();
}
std::vector<std::vector<TGraphErrors*> > RootItems::FillVecAvPH(std::vector<std::vector<TGraphErrors*> > graphVec){

    graphVec.resize(nRoc);
    for (uint8_t iRoc = 0; iRoc != nRoc; iRoc++) {
        for (uint8_t iMode = 0; iMode != 4; iMode++) {
            TGraphErrors * gr = new TGraphErrors;
            gr->SetName( Form("PulseHeightTime_ROC%i_NPix%i", iRoc, iMode) );
            gr->SetTitle( Form("Average Pulse Height ROC %i NPix %i", iRoc, iMode) );
            gr->GetXaxis()->SetTitle("Event Number");
            gr->GetYaxis()->SetTitle("Average Pulse Height (electrons)");
            gr->SetLineColor(HistColors[iMode]);
            gr->SetMarkerColor(HistColors[iMode]);
            gr->SetMinimum(0);
            gr->SetMaximum(70000);
            graphVec[iRoc].push_back(gr);
        }
    }
    return graphVec;
}
void RootItems::DrawSaveAvPH(uint8_t iroc){
    c1->cd();
    gAvgPH[iroc][0]->SetTitle(TString::Format("Average Pulse Height ROC%i", iroc));
    gAvgPH[iroc][0]->Draw("Ape");
    for (uint8_t iMode = 1; iMode != 4; iMode++)
        gAvgPH[iroc][iMode]->Draw("samepe");
    lPulseHeight->Draw("same");
    c1->SaveAs(OutDir + TString::Format("PulseHeightTime_ROC%i.gif", iroc));
}
void RootItems::AllocateArrAvPH(){
    dAvgPH2D = new double**[nRoc]; nAvgPH2D = new int**[nRoc];
    dAvgPH = new double*[nRoc]; nAvgPH = new int*[nRoc];
    for (uint8_t iRoc = 0; iRoc < nRoc; iRoc++){
        dAvgPH2D[iRoc] = new double*[PLTU::NCOL]; nAvgPH2D[iRoc] = new int*[PLTU::NCOL];
        dAvgPH[iRoc] = new double[4]; nAvgPH[iRoc] = new int[4];
        for (uint8_t iCol = 0; iCol < PLTU::NCOL; iCol++){
            dAvgPH2D[iRoc][iCol] = new double[PLTU::NROW]; nAvgPH2D[iRoc][iCol] = new int[PLTU::NROW];
        }
    }
}
vector<TH2F*> RootItems::FillVecResidual(vector<TH2F*> histVec, TString name, uint16_t xbin, float xmin, float xmax, uint16_t ybin, float ymin, float ymax){
    for (uint8_t iRoc = 0; iRoc != nRoc; iRoc++){
        TH2F * hist = new TH2F(Form(name, iRoc), Form(name, iRoc), xbin, xmin, xmax, ybin, ymin, ymax);
        histVec.push_back(hist);
    }
    return histVec;
}
void RootItems::DrawSaveResidual(uint8_t iroc, vector<TH2F*> histVec){
    c1->cd();
    gStyle->SetOptStat(1111);
    histVec[iroc]->Draw("colz");
    c1->SaveAs(OutDir+TString(histVec[iroc]->GetName()) + ".gif");
}
void RootItems::DrawSaveResidualProj(uint8_t iroc, vector<TH2F*> histVec, TString proj){
    c1->cd();
    if (proj == "X" or proj == "x"){
        histVec[iroc]->ProjectionX()->Draw();
        c1->SaveAs(OutDir+TString(hResidual[iroc]->GetName()) + "_X.gif");
    }
    else if (proj == "Y" or proj == "y"){
        histVec[iroc]->ProjectionY()->Draw();
        c1->SaveAs(OutDir+TString(hResidual[iroc]->GetName()) + "_Y.gif");
    }
}
/** Draw & Save */
void RootItems::DrawSaveOccupancy(uint8_t iroc, vector<TH2F*> histVec){

    c1->cd();
    histVec[iroc]->SetMinimum(0);
//    histVec[iroc].SetAxisRange(12,38,"X");
//    histVec[iroc].SetAxisRange(39,80,"Y");
    histVec[iroc]->Draw("colz");
    c1->SaveAs( OutDir+TString(histVec[iroc]->GetName()) + FileType);
}
void RootItems::DrawSaveOccupancy1DZ(uint8_t iroc){

    c1->cd();
    hOccupancy1DZ[iroc]->Draw("hist");
    if (hOccupancy1DZ[iroc]->GetEntries() > 0) c1->SetLogy(1);
    c1->SaveAs(OutDir+TString(hOccupancy1DZ[iroc]->GetName()) + FileType);
    c1->SetLogy(0);
}
void RootItems::DrawSaveOccupancyQuantile(uint8_t iroc){

    Double_t QProbability[1] = { 0.95 }; // Quantile positions in [0, 1]
//    Double_t QValue[1];
    hOccupancy1DZ[iroc]->GetQuantiles(1, QValue, QProbability); // saves threshold where 95% of the values are below in QValue
    if(QValue[0] > 1 && hOccupancy[iroc]->GetMaximum() > QValue[0])
        hOccupancy[iroc]->SetMaximum(QValue[0]);
    c1->cd();
    hOccupancy[iroc]->Draw("colz");
//    c1->SaveAs( OutDir+Form("Occupancy_ROC%i_Quantile.gif", iroc) );
    c1->SaveAs(OutDir + TString(hOccupancy[iroc]->GetName()) + "_Quantile" + FileType);

}
void RootItems::DrawSave3x3(uint8_t iroc){

    h3x3[iroc]->SetTitle( TString::Format("Occupancy Efficiency 3x3 ROC%i", iroc) );
    c1->cd();
    h3x3[iroc]->SetMinimum(0);
    h3x3[iroc]->SetMaximum(3);
    h3x3[iroc]->Draw("colz");
    c1->SaveAs(OutDir+TString(h3x3[iroc]->GetName()) + FileType);
}
void RootItems::DrawSave3x31DZ(uint8_t iroc){

    c1->cd();
    h3x31DZ[iroc]->Draw("hist");
    c1->SaveAs(OutDir+TString(h3x31DZ[iroc]->GetName()) + FileType);
}
void RootItems::DrawSaveAvPH2D(uint8_t iroc){

    c1->cd();
    hPulseHeightAvg2D[iroc]->SetMinimum(0);
    hPulseHeightAvg2D[iroc]->SetMaximum(100000);
    hPulseHeightAvg2D[iroc]->Draw("colz");
    c1->SaveAs(OutDir + hPulseHeightAvg2D[iroc]->GetName() + FileType);
}
void RootItems::DrawSaveTrackSlope(TH1F * slope){

    c1->cd();
    FitSlope(slope);
    slope->Draw();
    LegendSlope(slope);
    c1->SaveAs(OutDir + slope->GetName() + FileType);
}
