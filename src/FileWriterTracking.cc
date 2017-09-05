#include "FileWriterTracking.h"
#include "TInterpreter.h"

using std::cout; using std::string; using std::stringstream; using std::vector;

/** ============================
 CONSTRUCTOR
 =================================*/
FileWriterTracking::FileWriterTracking(string InFileName, uint8_t telescopeID, PSIFileReader * FR):
  nRoc(GetNumberOfROCS(telescopeID)), nHits(4), TelescopeID(telescopeID) {

  gInterpreter->GenerateDictionary("vector<vector<UShort_t> >");  // add vector<vector> to root dicts
  NewFileName = getFileName(InFileName);
  intree = ((PSIRootFileReader*) FR)->fTree;
  names = ((PSIRootFileReader*) FR)->fMacro;
  newfile = new TFile(NewFileName.c_str(), "RECREATE");
  br_dia_track_pos_x = new vector<float>;
  br_dia_track_pos_y = new vector<float>;
  br_dist_to_dia = new vector<float>;
  newtree = intree->CloneTree(0);
  br_charge_all.resize(nRoc);
  br_cluster_size = new vector<vector<uint16_t> >;
  br_cluster_size->resize(nRoc);
  br_clusters_per_plane.resize(nRoc);
  br_n_hits.resize(nRoc);
  br_track_x.resize(nRoc);
  br_track_y.resize(nRoc);
  br_smallest_hit_charge.resize(nRoc);
  br_smallest_hit_adc.resize(nRoc);
  br_smallest_hit_pos_col.resize(nRoc);
  br_smallest_hit_pos_row.resize(nRoc);
  br_residual_local_x.resize(nRoc);
  br_residual_local_y.resize(nRoc);
  br_residuals_x.resize(nRoc);
  br_residuals_y.resize(nRoc);
  addBranches();
}
FileWriterTracking::~FileWriterTracking() = default;

/** ============================
 AUXILIARY FUNCTIONS
 =================================*/
string FileWriterTracking::getFileName(string & InFileName){

  stringstream ss(InFileName);
  while (getline(ss, NewFileName, '/') != nullptr){}
  NewFileName.insert(unsigned(NewFileName.length() - 5), "_withTracks");
  return NewFileName;
}
void FileWriterTracking::addBranches(){

  newtree->Branch("hit_plane_bits", &br_hit_plane_bits);
  newtree->Branch("dia_track_x", &br_dia_track_pos_x);
  newtree->Branch("dia_track_y", &br_dia_track_pos_y);
  newtree->Branch("dist_to_dia", &br_dist_to_dia);
  newtree->Branch("chi2_tracks", &br_chi2);
  newtree->Branch("chi2_x", &br_chi2_x);
  newtree->Branch("chi2_y", &br_chi2_y);
  newtree->Branch("angle_x", &br_angle_x);
  newtree->Branch("angle_y", &br_angle_y);
  newtree->Branch("n_tracks", &br_n_tracks);
  newtree->Branch("n_clusters", &br_n_clusters);
  newtree->Branch("clusters_per_plane", &br_clusters_per_plane);
  newtree->Branch("cluster_plane", &br_cluster_plane);
  newtree->Branch("cluster_col", &br_cluster_col);
  newtree->Branch("cluster_row", &br_cluster_row);
  newtree->Branch("cluster_charge", &br_cluster_charge);
  newtree->Branch("cluster_xpos_tel", &br_cluster_xpos_tel);
  newtree->Branch("cluster_ypos_tel", &br_cluster_ypos_tel);
  newtree->Branch("cluster_xpos_local", &br_cluster_xpos_local);
  newtree->Branch("cluster_ypos_local", &br_cluster_ypos_local);
  newtree->Branch("n_hits", &br_n_hits);
  newtree->Branch("coincidence_map", &br_coincidence_map);
  newtree->Branch("residuals_x", &br_residuals_x);
  newtree->Branch("residuals_y", &br_residuals_y);
  newtree->Branch("cluster_size", &br_cluster_size);

  for (uint8_t iRoc = 0; iRoc != nRoc; iRoc++){
    TString branch_name_charge = TString::Format("charge_all_ROC%d", iRoc);
    newtree->Branch(branch_name_charge, &(br_charge_all[iRoc]));
    TString branch_name_residual_localX = TString::Format("residual_ROC%d_Local_X",iRoc);
    newtree->Branch(branch_name_residual_localX, &(br_residual_local_x[iRoc]));
    TString branch_name_residual_localY = TString::Format("residual_ROC%d_Local_Y",iRoc);
    newtree->Branch(branch_name_residual_localY, &(br_residual_local_y[iRoc]));
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
  }
}

void FileWriterTracking::fillTree(){
  newtree->Fill();
}

void FileWriterTracking::saveTree(){

  newfile->cd();
  newtree->Write();
  if (names != nullptr)
    names->Write();
  newfile->Write();
  newfile->Close();
  delete newfile;
}

void FileWriterTracking::clearVectors(){

  br_cluster_plane.clear();
  br_cluster_col.clear();
  br_cluster_row.clear();
  br_cluster_charge.clear();
  br_cluster_xpos_tel.clear();
  br_cluster_ypos_tel.clear();
  br_cluster_xpos_local.clear();
  br_cluster_ypos_local.clear();
  for (uint8_t iRoc = 0; iRoc != nRoc; iRoc++) {
    br_charge_all[iRoc]->clear();
    br_cluster_size->at(iRoc).clear();
    br_residual_local_x[iRoc]->clear();
    br_residual_local_y[iRoc]->clear();
    br_track_x[iRoc]->clear();
    br_track_y[iRoc]->clear();
    br_smallest_hit_charge[iRoc]->clear();
    br_smallest_hit_adc[iRoc]->clear();
    br_smallest_hit_pos_col[iRoc]->clear();
    br_smallest_hit_pos_row[iRoc]->clear();
  }
  br_dia_track_pos_x->clear();
  br_dia_track_pos_y->clear();
  br_dist_to_dia->clear();
}

