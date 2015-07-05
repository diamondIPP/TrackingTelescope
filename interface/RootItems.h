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

/** ============================
 ROOTITEMS CLASS
 =================================*/

class RootItems {

/** declare all the histograms, graphs, etc. */
private:

    /** constants */
    uint16_t const nRoc;
    TString const PlotsDir;
    TString const OutDir;

    /** canvases */
    TCanvas * c2;

    /** tracking */
    TH1F * hTrackSlopeX;
    TH1F * hTrackSlopeY;
    TF1 * fGauss;
    TLegend * lFitGauss;

    /** occupancy */
    std::vector<TH2F*> hOccupancy;
    std::vector<TH2F*> hOccupancyLowPH;
    std::vector<TH2F*> hOccupancyHighPH;

    /** cluster hits */
    std::vector<TH1F*> hNHitsPerCluster;
    std::vector<TH1F*> hNClusters;

    /** coincidence map */
    TH1F * hCoincidenceMap;

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

    /** ============================
     AUXILIARY FUNCTIONS
     =================================*/
    void FitSlope(TH1F * histo);
    void LegendSlope(TH1F * histo);
    std::vector<TH2F*> FillVectorTH2F(std::vector<TH2F*> histo, const char * name);
    std::vector<TH1F*> FillVectorTH1F(std::vector<TH1F*> histo, const char * name);
    void DrawSaveTH1F(std::vector<TH1F*> histo, uint8_t iroc, TCanvas & c, const char * xTit, const char * yTit);
    void PrepCoincidenceHisto();
    void DrawSaveCoincidence();



 };




#endif // RootItems_h
