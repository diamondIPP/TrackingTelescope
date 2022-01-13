#include "PSIRootFileReader.h"

#include <iostream>
#include <string>
#include <utility>
#include <cstdint>

using namespace std;

PSIRootFileReader::PSIRootFileReader(string in_file_name, bool const only_align, bool track_only_telescope):
  PSIFileReader(track_only_telescope), fFileName(move(in_file_name)), fOnlyAlign(only_align) {
    if (!OpenFile()) {
        std::cerr << "ERROR: cannot open input file: " << fFileName << std::endl;
    throw;
    }
}


PSIRootFileReader::~PSIRootFileReader ()
{
  Clear();
  CloseFile();
  //  delete fTree;
  // Crashed when uncommented. Live with the memleak for now
  //delete fRootFile;
}


bool PSIRootFileReader::OpenFile ()
{
    f_n_hits = 0;
    cout << "Open File " << fFileName << endl;
    fRootFile = new TFile(fFileName.c_str(), "READ");

    if (!fRootFile->IsOpen()) { return false; }

    fTree = dynamic_cast<TTree*>(fRootFile->Get("tree"));
    if(fRootFile->FindKey("region_information") != nullptr) {
        fMacro = dynamic_cast<TMacro *>(fRootFile->Get("region_information"));
    }

    fAtEntry = 0;
    fNEntries = int(fTree->GetEntries());

    if (fNEntries <= 0) return false;

    // Set Branch Addresses
    fTree->SetBranchAddress("n_hits_tot", &f_n_hits);
    fTree->SetBranchAddress("event_number", &f_event_number);
    fTree->SetBranchAddress("time", &f_time);

    fTree->SetBranchAddress("plane", f_plane);
    fTree->SetBranchAddress("col", f_col);
    fTree->SetBranchAddress("row", f_row);
    fTree->SetBranchAddress("adc", f_adc);
    fTree->SetBranchAddress("charge", f_charge);
    if (fTree->FindBranch(GetSignalBranchName() ))
        fTree->SetBranchAddress(GetSignalBranchName(), f_signal);
    return true;
}

void PSIRootFileReader::CloseFile() {
    if(fRootFile->IsOpen()) {
        fRootFile->Close();
        cout << "Close File" << endl;
    } else
        cout << "File is closed" << endl;
}

void PSIRootFileReader::ResetFile ()
{
    CloseFile();
    cout << "Reset File" << endl;
    OpenFile();
}

int PSIRootFileReader::GetNextEvent ()
{
    Clear();
    if ( !fOnlyAlign ){
        if (fTree->GetBranch(GetSignalBranchName() )){
            AddSignal(f_signal, f_n_hits);
        }
    }

    for (int i = 0; i != fNPlanes; ++i) {
        fPlaneMap.emplace(i, PLTPlane());
        fPlaneMap[i].SetROC(i);
    }

    if (fAtEntry == fNEntries) {
        return -1;
    }

    fTree->GetEntry(fAtEntry);

    fAtEntry++;
    if (f_n_hits > 255) { cout << endl<< "f_plane->size() = " << f_n_hits << endl; }

    for (auto i_hit = 0; i_hit != f_n_hits; i_hit++){
        uint8_t roc = f_plane[i_hit];
        uint8_t col = f_col[i_hit];
        uint8_t row = f_row[i_hit];
        int16_t adc = f_adc[i_hit];

        if (!IsPixelMasked( 1*100000 + roc*10000 + col*100 + row)){
            auto * Hit = new PLTHit(1, roc, col, row, adc);

            /** Gain calibration */
            fGainCal.SetCharge(*Hit);
            f_charge[i_hit] = Hit->Charge();  // overwrite empty charge values...

            /** Alignment */
            fAlignment.AlignHit(*Hit);
            fHits.push_back(Hit);
            fPlaneMap[Hit->ROC()].AddHit(Hit);
            if ( fOnlyAlign ) {
                for (uint8_t i = 0; i !=roc+1; i++) {
                    if (fPlaneMap[i].NHits() == 0) return 0;
                }
            } // CHECKS THAT THERE WERE HITS IN THE PREVIOUS ROCS IF NOT RETURN 0
        }
    }

    /** Loop over all planes and clusterize each one, then add each plane to the correct telescope (by channel number) */
    for (auto & it : fPlaneMap){
        it.second.Clusterize(PLTPlane::kClustering_AllTouching, PLTPlane::kFiducialRegion_All);
        AddPlane( &(it.second) );
    }

    /** If we are doing single plane-efficiencies:
        Just send all events to the tracking and sort it out there */
    if (DoingSinglePlaneEfficiency()) { // TRUE FOR 1ST LOOP ALIGNMENT
        RunTracking( *((PLTTelescope*) this));
    }

    /** Otherwise require exactly one cluster per plane */
    else{
        switch (fTrackingAlgorithm) {
            case kTrackingAlgorithm_ETH:
                if (HaveOneCluster(4) && NClusters() != NPlanes() && (HitPlaneBits() & 15) == pow(2, 4) - 1){
//                    cout << "Event has the required conditions (ETH tracking): NClusters: " << NClusters() << " and HitPlaneBits (15): " << (HitPlaneBits() & 15) << endl;
                    RunTracking(*((PLTTelescope*)this));
                }
            case kTrackingAlgorithm_6PlanesHit:
                if (NClusters() == NPlanes() && HitPlaneBits() == pow(2, NPlanes() ) - 1){
//                    cout << "Event has the required conditions (All planes for tracking): NClusters: " << NClusters() << " and HitPlaneBits (127): " << HitPlaneBits() << endl;
                    RunTracking( *((PLTTelescope*) this));
                }
//            default:
//                cout << "Entered the default for tracking " << endl;
//                if (NClusters() == NPlanes() && HitPlaneBits() == pow(2, NPlanes() ) - 1){
//                    cout << "Event has the required conditions: NClusters: " << NClusters() << " and HitPlaneBits (127): " << HitPlaneBits() << endl;
//                    RunTracking( *((PLTTelescope*) this));
//                }
            case kTrackingAlgorithm_NoTracking:break;
            case kTrackingAlgorithm_01to2_All:break;
            case kTrackingAlgorithm_2PlaneTracks_All:break;
        }
    }
    return 0;

}