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

/** ============================
 ROOTITEMS CLASS
 =================================*/

class RootItems {

/** declare all the histograms, graphs, etc. */
private:

    /** tracking */
    TH1F * hTrackSlopeX;
    TH1F * hTrackSlopeY;
    TF1 * fGauss;
    TLegend * lFitGauss;

public:

    /** ============================
     CONSTRUCTOR
     =================================*/
    RootItems();
    ~RootItems();


    /** ============================
     GET-FUNCTIONS
     =================================*/
    TH1F * TrackSlopeX();
    TH1F * TrackSlopeY();
    TLegend * FitGauss();


    /** ============================
     HELPER FUNCTIONS
     =================================*/
    void FitSlope(TH1F * histo);
    void LegendSlope(TH1F * histo);



 };








#endif // RootItems_h
