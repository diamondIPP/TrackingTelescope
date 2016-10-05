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
    newfile = new TFile(NewFileName.c_str(), "RECREATE");
    newtree = intree->CloneTree(0);
//    br_pulse_heights_all.resize(nRoc);
//    for(uint8_t iRoc = 0; iRoc != nRoc; iRoc++){
//        br_pulse_heights_all[iRoc].resize(nHits);
//    }
    br_charge_all.resize(nRoc);
    br_cluster_size.resize(nRoc);
    br_clusters_per_plane.resize(nRoc);
    br_cluster_pos_telescope_x.resize(nRoc);
    br_cluster_pos_telescope_y.resize(nRoc);
    br_cluster_pos_local_x.resize(nRoc);
    br_cluster_pos_local_y.resize(nRoc);
    br_cluster_col.resize(nRoc);
    br_cluster_row.resize(nRoc);
    br_track_x.resize(nRoc);
    br_track_y.resize(nRoc);
    br_smallest_hit_charge.resize(nRoc);
    br_smallest_hit_adc.resize(nRoc);
    br_smallest_hit_pos_col.resize(nRoc);
    br_smallest_hit_pos_row.resize(nRoc);
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
    newtree->Branch("angle_x", &br_angle_x);
    newtree->Branch("angle_y", &br_angle_y);
    newtree->Branch("n_tracks", &br_n_tracks);
    newtree->Branch("n_clusters", &br_n_clusters);
    newtree->Branch("clusters_per_plane", &br_clusters_per_plane);
    // new branches: DA
    newtree->Branch("coincidence_map", &br_coincidence_map);

    for (uint8_t iRoc = 0; iRoc != nRoc; iRoc++){
        TString branch_name_charge = TString::Format("charge_all_ROC%d", iRoc);
        newtree->Branch(branch_name_charge, &(br_charge_all[iRoc]));
        TString branch_name_cluster_size = TString::Format("cluster_size_ROC%d", iRoc);
        newtree->Branch(branch_name_cluster_size, &(br_cluster_size[iRoc]));
        TString branch_name_cluster_telescopeX = TString::Format("cluster_pos_ROC%d_Telescope_X", iRoc);
        newtree->Branch(branch_name_cluster_telescopeX, &(br_cluster_pos_telescope_x[iRoc]));
        TString branch_name_cluster_telescopeY = TString::Format("cluster_pos_ROC%d_Telescope_Y", iRoc);
        newtree->Branch(branch_name_cluster_telescopeY, &(br_cluster_pos_telescope_y[iRoc]));
        TString branch_name_cluster_localX = TString::Format("cluster_pos_ROC%d_Local_X",iRoc);
        newtree->Branch(branch_name_cluster_localX, &(br_cluster_pos_local_x[iRoc]));
        TString branch_name_cluster_localY = TString::Format("cluster_pos_ROC%d_Local_Y",iRoc);
        newtree->Branch(branch_name_cluster_localY, &(br_cluster_pos_local_y[iRoc]));
        TString branch_name_cluster_row = TString::Format("cluster_row_ROC%d",iRoc);
        newtree->Branch(branch_name_cluster_row, &(br_cluster_row[iRoc]));
        TString branch_name_cluster_col = TString::Format("cluster_col_ROC%d",iRoc);
        newtree->Branch(branch_name_cluster_col, &(br_cluster_col[iRoc]));
        TString branch_name_track_x = TString::Format("track_x_ROC%d",iRoc);
        newtree->Branch(branch_name_track_x, &(br_track_x[iRoc]));
        TString branch_name_track_y = TString::Format("track_y_ROC%d",iRoc);
        newtree->Branch(branch_name_track_y, &(br_track_y[iRoc]));
        TString branch_name_smallest_charge = TString::Format("smallest_clust_hit_charge_ROC%d",iRoc);
        newtree->Branch(branch_name_smallest_charge, &(br_smallest_hit_charge[iRoc]));
        TString branch_name_smallest_adc = TString::Format("smallest_clust_hit_adc_ROC%d",iRoc);
        newtree->Branch(branch_name_smallest_adc, &(br_smallest_hit_adc[iRoc]));
        TString branch_name_smallest_col = TString::Format("smallest_clust_hit_col_ROC%d",iRoc);
        newtree->Branch(branch_name_smallest_col, &(br_smallest_hit_pos_col[iRoc]));
        TString branch_name_smallest_row = TString::Format("smallest_clust_hit_row_ROC%d",iRoc);
        newtree->Branch(branch_name_smallest_row, &(br_smallest_hit_pos_row[iRoc]));

//        for(size_t iHits = 1; iHits < nHits; iHits++){
//            TString branch_name_RocPulseHeights = TString::Format("pulse_height_ROC%d_%d_cluster",iRoc,iHits);
//            newtree->Branch(branch_name_RocPulseHeights,&(br_pulse_heights_all[iRoc][iHits-1]));
//        }
//        TString branch_name_RocPulseHeightsLast = TString::Format("pulse_height_ROC%d_More%d_cluster",iRoc,nHits);
//        newtree->Branch(branch_name_RocPulseHeightsLast,&(br_pulse_heights_all[iRoc][nHits-1]));
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
        br_cluster_size[iRoc]->clear();
        br_cluster_pos_telescope_x[iRoc]->clear();
        br_cluster_pos_telescope_y[iRoc]->clear();
        br_cluster_pos_local_x[iRoc]->clear();
        br_cluster_pos_local_y[iRoc]->clear();
        br_cluster_col[iRoc]->clear();
        br_cluster_row[iRoc]->clear();
        br_track_x[iRoc]->clear();
        br_track_y[iRoc]->clear();
        br_smallest_hit_charge[iRoc]->clear();
        br_smallest_hit_adc[iRoc]->clear();
        br_smallest_hit_pos_col[iRoc]->clear();
        br_smallest_hit_pos_row[iRoc]->clear();
//        for (uint8_t iHits = 0; iHits != nHits; iHits++){
//            br_pulse_heights_all[iRoc][iHits]->clear();
//        }
    }
}

