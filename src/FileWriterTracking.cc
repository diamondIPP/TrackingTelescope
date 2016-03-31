#include "FileWriterTracking.h"
#include "GetNames.h"

using namespace std;

/** ============================
 CONSTRUCTOR
 =================================*/
FileWriterTracking::FileWriterTracking(string InFileName, uint8_t telescopeID, PSIFileReader * FR):
    nRoc(GetNumberOfROCS(telescopeID)), nHits(4) {

    NewFileName = getFileName(InFileName);
    intree = ((PSIRootFileReader*) FR)->fTree;
    names = ((PSIRootFileReader*) FR)->fMacro;
    if(names){
        names->Print();
    }
    newfile = new TFile(NewFileName.c_str(), "RECREATE");
    newtree = intree->CloneTree(0);
    br_pulse_heights_all.resize(nRoc);
    for(uint8_t iRoc = 0; iRoc != nRoc; iRoc++){
        br_pulse_heights_all[iRoc].resize(nHits);
    }
    br_charge_all.resize(nRoc);
    br_clusters_per_plane.resize(nRoc);
    br_cluster_pos_x.resize(nRoc);
    br_cluster_pos_y.resize(nRoc);
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

    for (uint8_t iRoc = 0; iRoc != nRoc; iRoc++){
        TString branch_name_charge = TString::Format("charge_all_ROC%d", iRoc);
        newtree->Branch(branch_name_charge, &(br_charge_all[iRoc]));
        TString branch_name_clusterX = TString::Format("cluster_pos_ROC%d_X", iRoc);
        newtree->Branch(branch_name_clusterX, &(br_cluster_pos_x[iRoc]));
        TString branch_name_clusterY = TString::Format("cluster_pos_ROC%d_Y", iRoc);
        newtree->Branch(branch_name_clusterY, &(br_cluster_pos_y[iRoc]));
        for(size_t iHits = 1; iHits < nHits; iHits++){
            TString branch_name_RocPulseHeights = TString::Format("pulse_height_ROC%d_%d_cluster",iRoc,iHits);
            newtree->Branch(branch_name_RocPulseHeights,&(br_pulse_heights_all[iRoc][iHits-1]));
        }
        TString branch_name_RocPulseHeightsLast = TString::Format("pulse_height_ROC%d_More%d_cluster",iRoc,nHits);
        newtree->Branch(branch_name_RocPulseHeightsLast,&(br_pulse_heights_all[iRoc][nHits-1]));
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
        br_cluster_pos_x[iRoc]->clear();
        br_cluster_pos_y[iRoc]->clear();
        for (uint8_t iHits = 0; iHits != nHits; iHits++){
            br_pulse_heights_all[iRoc][iHits]->clear();
        }
    }
}

