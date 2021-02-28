#include "TestPlaneEfficiencySilicon.h"

using namespace std;


/** TestPlaneEfficiencySilicon
  o) Consider one plane to be the plane under test
  o) Require at least one hit in both silicon planes
  o) Check if we have a hit in the other planes */

int TestPlaneEfficiencySilicon (std::string const InFileName, TFile * out_f,
                                 TString const RunNumber, int telescopeID)
{
    // TODO: Currently hard coded number of planes
    // Was fine as we only did pixel-analysis (where this is called)
    // with a six-plane telescope before.
    // Fix this!

    if (DEBUG) cout << "DEBUG: Entering TestPlaneEfficiencySilicon" << endl;

    gStyle->SetOptStat(0);
    TString const PlotsDir = "plots/";
    TString const OutDir = PlotsDir + RunNumber + "/";


    if (DEBUG) cout << "DEBUG: Initializing BinaryFileReader" << endl;

    /** Initialize Reader */
    PSIFileReader * FR;

    if (IsROOTFile(InFileName)){
        FR = new PSIRootFileReader(InFileName, 0, false);
    }
    else {
        FR = new PSIBinaryFileReader(InFileName);
        ((PSIBinaryFileReader*) FR)->CalculateLevels(OutDir);
    }

    FR->GetAlignment()->SetErrors(telescopeID);

    /** Apply Masking */
    // TODO: More dynamic selection of masking for this
    FR->ReadPixelMask("outerPixelMask_forSiEff.txt");

    /** numerators and denominators for efficiency calculation */
    vector<int> nums(6);
    vector<int> denoms(6);
    for (int i = 0; i != 6; i++){
        nums[i]   = 0;
        denoms[i] = 0;
    }

    int n_events = 0;

    /** Event Loop */
    for (int ievent = 0; FR->GetNextEvent() >= 0; ++ievent) {

        n_events++;

        /** print progress */
        if (ievent % 10000 == 0) {
          cout << "Processing event: " << ievent << endl;
        }

        /** Initializes all planes as un-hit */
        vector<bool> plane_is_hit;
        for (int i=0;i!=6;i++)
            plane_is_hit.push_back(false);

        /** then loop and see where we have a hit */
        for (uint32_t ihit = 0; ihit != FR->NHits(); ++ihit)
            plane_is_hit[ FR->Hit(ihit)->ROC() ] = true;

        /** Check for Coincidence of silicon hits */
        if (plane_is_hit[0] && plane_is_hit[5]){

            /** Increment denominators */
            for (int i = 0; i != 6; i++)
                denoms[i]++;

            /** Increment numerators */
            for (int i = 0; i != 6; i++)
                if (plane_is_hit[i])
                    nums[i]++;
        }

    } // End of Event Loop

    out_f->cd();

    /** Store the numerators and denominators to the file */
    for (int i=0; i!=6; i++){

        /** numerator */
        TParameter<int> n;
        n.SetVal( nums[i]);
        n.Write( Form("SiliconEfficiencyNumeratorROC%i",i));

        /** denominator */
        TParameter<int> d;
        d.SetVal( denoms[i]);
        d.Write( Form("SiliconEfficiencyDenominatorROC%i",i));

    }

    if (DEBUG) cout << "DEBUG: Leaving TestPlaneEfficiencySilicon" << endl;

    delete FR;

    return n_events;
}


