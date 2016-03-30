#include "FileWriterTracking.h"
#include "GetNames.h"

using namespace std;

/** ============================
 CONSTRUCTOR
 =================================*/
FileWriterTracking::FileWriterTracking(string InFileName, uint8_t telescopeID, PSIFileReader * FR):
    nRoc(GetNumberOfROCS(telescopeID)) {

    NewFileName = getFileName(InFileName);
    intree = ((PSIRootFileReader*) FR)->fTree;
    names = ((PSIRootFileReader*) FR)->fMacro;
    if(names){
        names->Print();
    }
    newfile = new TFile(NewFileName.c_str(), "RECREATE");
    newtree = intree->CloneTree(0);
    // Start vector if declared as a pointer: DA
    //br_pulse_height = new std::vector<float>;
    br_charge_all.resize(nRoc);
    br_clusters_per_plane.resize(nRoc);
//    br_cluster_pos_x.resize(nRoc);
//    br_cluster_pos_y.resize(nRoc);
    addBranches();

}
FileWriterTracking::~FileWriterTracking(){

}

/** ============================
 AUXILIARY FUNCTIONS
 =================================*/
string FileWriterTracking::getFileName(string InFileName){
    stringstream ss(InFileName);
    while (getline(ss, NewFileName, '/')){}
    NewFileName.insert(int(NewFileName.length()-5), "_withTracks");
    return NewFileName;
}
void FileWriterTracking::addBranches(){
    newtree->Branch("hit_plane_bits", &br_hit_plane_bits);
    newtree->Branch("diam1_track_x", &br_diam1_track_x);
    newtree->Branch("diam1_track_y", &br_diam1_track_y);
    newtree->Branch("diam2_track_x", &br_diam2_track_x);
    newtree->Branch("diam2_track_y", &br_diam2_track_y);
    newtree->Branch("dist_to_dia1", &br_dist_to_dia1);
    newtree->Branch("dist_to_dia2", &br_dist_to_dia2);
    newtree->Branch("chi2_tracks", &br_chi2);
    newtree->Branch("chi2_x", &br_chi2_x);
    newtree->Branch("chi2_y", &br_chi2_y);
    newtree->Branch("slope_x", &br_slope_x);
    newtree->Branch("slope_y", &br_slope_y);
    newtree->Branch("n_tracks", &br_n_tracks);
    newtree->Branch("n_clusters", &br_n_clusters);
    newtree->Branch("clusters_per_plane", &br_clusters_per_plane);
//    newtree->Branch("cluster_pos_x", &br_cluster_pos_x);
//    newtree->Branch("cluster_pos_y", &br_cluster_pos_y);
//    newtree->Branch("test", &br_test);
    // Add branch pulse height to new tree: DA
    //newtree->Branch("pulse_height",&br_pulse_height);
    for (uint8_t iRoc = 0; iRoc != nRoc; iRoc++){
        TString branch_name = TString::Format("charge_all_ROC%d", iRoc);
        newtree->Branch(branch_name, &(br_charge_all[iRoc]));
    }
}
void FileWriterTracking::fillTree(){
    newtree->Fill();
}
void FileWriterTracking::saveTree(){
    newfile->cd();
    newtree->Write();
    if (names) names->Write();
    newfile->Write();
    newfile->Close();
    delete newfile;

}
void FileWriterTracking::clearVectors(){
    for (uint8_t iRoc = 0; iRoc !=nRoc; iRoc++){
        br_charge_all[iRoc]->clear();
//        br_cluster_pos_x[iRoc].clear();
//        br_cluster_pos_y[iRoc].clear();
    }
    // Clear vector of pulse height: DA
    //br_pulse_height.clear();
}

