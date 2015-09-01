#ifndef RootItems_h
#define RootItems_h

#include <iostream>
#include <vector>
#include <string>
#include <ctime>

#include "TLegend.h"
#include "TLegendEntry.h"
#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH3F.h"
#include "TProfile2D.h"
#include "TParameter.h"
#include "TTree.h"
#include "TText.h"
#include "TPad.h"
#include "TF2.h"
#include "TPaveStats.h"

#include "PSIRootFileReader.h"
#include "PSIBinaryFileReader.h"
#include "GetNames.h"
#include "PLTU.h"

/** ============================
 ROOTITEMS CLASS
 =================================*/

class RootItems {

/** declare all the histograms, graphs, etc. */
private:

    /** constants */
    uint16_t const nRoc, nSig;
    TString const PlotsDir;
    TString const OutDir;
    TString const FileType;
    const uint8_t HistColors[4];
    const uint8_t maxChi2;

    /** canvases */
    TCanvas * c1, * c2;

    /** tracking */
    TH1F * hTrackSlopeX;
    TH1F * hTrackSlopeY;
    TF1 * fGauss;
    TLegend * lFitGauss;

    /** occupancy */
    std::vector<TH2F*> hOccupancy;
    std::vector<TH1F*> hOccupancy1DZ;
    Double_t QValue[1];
    std::vector<TH2F*> h3x3;
    std::vector<TH1F*> h3x31DZ;
    std::vector<TH2F*> hOccupancyLowPH;
    std::vector<TH2F*> hOccupancyHighPH;

    /** cluster hits */
    std::vector<TH1F*> hNHitsPerCluster;
    std::vector<TH1F*> hNClusters;

    /** pulse height */
    std::vector<std::vector<TH1F*> > hPulseHeight;
    std::vector<std::vector<TH1F*> > hPulseHeightLong;
    std::vector<std::vector<TH1F*> > hPulseHeightOffline;
    TLegend * lPulseHeight;
    TLegend * lPHMean;
    TLegend * lRatio;
    std::vector<std::vector<TGraphErrors*> > gAvgPH;
    double *** dAvgPH2D;
    int *** nAvgPH2D;
    double ** dAvgPH;
    int ** nAvgPH;
    std::vector<TH2F*> hPulseHeightAvg2D;

    /** coincidence map */
    TH1F * hCoincidenceMap;

    /** chi2 */
    TH1F * hChi2;
    TH1F * hChi2X;
    TH1F * hChi2Y;

    /** residuals */
    std::vector<TH2F*> hResidual;
    std::vector<TH2F*> hResidualXdY;
    std::vector<TH2F*> hResidualYdX;

    /** signal distribution */
    std::vector<TProfile2D*> hSignalDistribution;


public:

    /** ============================
     CONSTRUCTOR
     =================================*/
    RootItems(uint8_t telescopeID, TString const RunNumber);
    ~RootItems();


    /** ============================
     GET-FUNCTIONS
     =================================*/
    TH1F * TrackSlopeX() { return hTrackSlopeX; }
    TH1F * TrackSlopeY() { return hTrackSlopeY; }
    TH1F * CoincidenceMap() { return hCoincidenceMap; }
    TLegend * FitGauss() { return lFitGauss; }
    std::vector<TH2F*> Occupancy() { return hOccupancy; }
    std::vector<TH2F*> OccupancyLowPH() { return hOccupancyLowPH; }
    std::vector<TH2F*> OccupancyHighPH() { return hOccupancyHighPH; }
    std::vector<TH1F*> nHitsPerCluster() { return hNHitsPerCluster; }
    std::vector<TH1F*> nClusters() { return hNClusters; }
    std::vector<std::vector<TH1F*> > PulseHeight() { return hPulseHeight; }
    std::vector<std::vector<TH1F*> > PulseHeightLong() { return hPulseHeightLong; }
    std::vector<std::vector<TH1F*> > PulseHeightOffline() { return hPulseHeightOffline; }
    std::vector<TProfile2D*> SignalDisto() { return hSignalDistribution; }
    TLegend * legPH() { return lPulseHeight; }
    TLegend * legPHMean() { return lPHMean; }
    TH1F * Chi2() { return hChi2; }
    TH1F * Chi2X() { return hChi2X; }
    TH1F * Chi2Y() { return hChi2Y; }
    std::vector<std::vector<TGraphErrors*> > AveragePH() { return gAvgPH; }
    double *** dAveragePH2D() { return dAvgPH2D; }
    int *** nAveragePH2D() { return nAvgPH2D; }
    double ** dAveragePH() { return dAvgPH; }
    int ** nAveragePH() { return nAvgPH; }
    std::vector<TH2F*> Residual() { return hResidual; }
    std::vector<TH2F*> ResidualXdY() { return hResidualXdY; }
    std::vector<TH2F*> ResidualYdX() { return hResidualYdX; }
    std::vector<TH1F*> Occupancy1DZ() {return hOccupancy1DZ; }
    Double_t * QuantileValue() { return QValue; }
    std::vector<TH2F*> Eff3x3() { return h3x3; }
    std::vector<TH2F*> PulseHeightAv2D() { return hPulseHeightAvg2D; }
    uint8_t NRoc() { return nRoc; }
    uint8_t NSig() { return nSig; }
    TString getOutDir() { return OutDir; }
    TString getPlotsDir() { return PlotsDir; }


    /** ============================
     SET-FUNCTIONS
     =================================*/
     void setOccupancy1DZ(TH1F * histo, uint8_t iroc) {hOccupancy1DZ[iroc] = histo; }
     void set3x3(TH2F * histo, uint8_t iroc) {h3x3[iroc] = histo; }
     void set3x31DZ(TH1F * histo, uint8_t iroc) {h3x31DZ[iroc] = histo; }


     /** ============================
     MAIN FUNCTIONS
     =================================*/
     void SaveAllHistos();


    /** ============================
     AUXILIARY FUNCTIONS
     =================================*/
    /** Fill vectors */
    std::vector<TProfile2D*> FillSignalDisto();
    void FitSlope(TH1F * histo);
    TH1F * FormatSlopeHisto(TString, uint16_t, float);
    void LegendSlope(TH1F * histo, TString);
    std::vector<TH2F*> FillVectorTH2F(std::vector<TH2F*> histo, const char * name);
    std::vector<TH1F*> FillVectorTH1F(std::vector<TH1F*> histo, const char * name);
    std::vector<std::vector<TH1F*> > FillVectorPH(std::vector<std::vector<TH1F*> >, TString name, uint32_t maxPH);
    void FillAvPH2D(uint8_t);
    void PrepCoincidenceHisto();
    void FormatPHHisto(std::vector<std::vector<TH1F*> >);
    void FormatLegendsPH();
    void FillLegendsPH(uint8_t iroc, std::vector<std::vector<TH1F*> > histVec);
    void ClearLegendsPH();
    void DrawSaveChi2(TH1F*, TString);
    std::vector<std::vector<TGraphErrors*> > FillVecAvPH(std::vector<std::vector<TGraphErrors*> >);
    void AllocateArrAvPH();
    std::vector<TH2F*> FillVecResidual(std::vector<TH2F*>, TString name, uint16_t, float, float, uint16_t, float, float);
    /** Draw & Save */
    void DrawSaveCoincidence();
    void DrawSavePH(uint8_t iroc, std::vector<std::vector<TH1F*> > histVec, TString title, TString saveName);
    void DrawSaveAvPH(uint8_t);
    void DrawSaveTH1F(std::vector<TH1F*> histo, uint8_t iroc, const char * xTit, const char * yTit);
    void DrawSaveResidual(uint8_t, vector<TH2F*>);
    void DrawSaveResidualProj(uint8_t, vector<TH2F*>, TString);
    void DrawSaveOccupancy(uint8_t, vector<TH2F*>);
    void DrawSaveOccupancy1DZ(uint8_t);
    void DrawSaveOccupancyQuantile(uint8_t);
    void DrawSave3x3(uint8_t);
    void DrawSave3x31DZ(uint8_t);
    void DrawSaveAvPH2D(uint8_t);
    void DrawSaveTrackSlope(TH1F*);
    void DrawSaveSignalDisto();
 };




#endif // RootItems_h
