#include "FileWriterTracking.h"
#include "Utils.h"

using std::cout; using std::string; using std::stringstream; using std::vector; using std::endl;

/** ============================
 CONSTRUCTOR
 =================================*/
FileWriterTracking::FileWriterTracking(string InFileName, PSIFileReader * FR):
  n_rocs_(GetNPlanes()), n_duts_(GetNDUTs()), FR_(FR) {

  NewFileName = getFileName(InFileName);
  intree = ((PSIRootFileReader*) FR)->fTree;
  names = ((PSIRootFileReader*) FR)->fMacro;
  newfile = new TFile(NewFileName.c_str(), "RECREATE");
  newtree = intree->CloneTree(0);

  /** init arrays */
  br_dia_track_pos_x = new float[n_duts_];
  br_dia_track_pos_y = new float[n_duts_];
  br_dia_track_pos_x_loc = new float[n_duts_];
  br_dia_track_pos_y_loc = new float[n_duts_];
  br_dist_to_dia = new float[n_duts_];
  br_n_hits = new uint16_t[n_rocs_];
  br_n_clusters = new uint8_t[n_rocs_];
  br_sres = new float[n_rocs_];
  br_sres_x = new float[n_rocs_];
  br_sres_y = new float[n_rocs_];

  /** init vectors */
  br_residuals_x = new vector<vector<float> >;
  br_residuals_y = new vector<vector<float> >;
  br_residuals = new vector<vector<float> >;
  br_track_x = new vector<vector<float> >;
  br_track_y = new vector<vector<float> >;

  br_cluster_size = new vector<vector<uint16_t> >;
  br_cluster_col = new vector<vector<uint16_t> >;
  br_cluster_row = new vector<vector<uint16_t> >;
  br_cluster_xpos_tel = new vector<vector<float> >;
  br_cluster_ypos_tel = new vector<vector<float> >;
  br_cluster_xpos_local = new vector<vector<float> >;
  br_cluster_ypos_local = new vector<vector<float> >;

  br_cluster_charge = new vector<vector<float> >;
  br_aligned = new vector<bool>;
  is_aligned = new vector<bool>;
  resizeVectors();
  addBranches();
}

/** ============================
 AUXILIARY FUNCTIONS
 =================================*/
string FileWriterTracking::getFileName(string & InFileName){

  stringstream ss(InFileName);
  while (getline(ss, NewFileName, '/') ){}
  NewFileName.insert(unsigned(NewFileName.length() - 5), "_withTracks");
  return NewFileName;
}
void FileWriterTracking::addBranches(){

  newtree->Branch("hit_plane_bits", &br_hit_plane_bits);
  /** arrays */
  newtree->Branch("dia_track_x", br_dia_track_pos_x, Form("dia_track_x[%d]/F", n_duts_));
  newtree->Branch("dia_track_y", br_dia_track_pos_y, Form("dia_track_y[%d]/F", n_duts_));
  newtree->Branch("dia_track_x_local", br_dia_track_pos_x_loc, Form("dia_track_x_local[%d]/F", n_duts_));
  newtree->Branch("dia_track_y_local", br_dia_track_pos_y_loc, Form("dia_track_Y_local[%d]/F", n_duts_));
  newtree->Branch("dist_to_dia", br_dist_to_dia, Form("dist_to_dia[%d]/F", n_duts_));
  newtree->Branch("n_hits", br_n_hits, Form("n_hits[%d]/s", n_rocs_));
  newtree->Branch("n_clusters", br_n_clusters, Form("n_clusters[%d]/b", n_rocs_));
  newtree->Branch("sres", br_sres, Form("sres[%d]/F", n_rocs_));
  newtree->Branch("sres_x", br_sres_x, Form("sres_x[%d]/F", n_rocs_));
  newtree->Branch("sres_y", br_sres_y, Form("sres_y[%d]/F", n_rocs_));
  /** vectors */
  newtree->Branch("chi2_tracks", &br_chi2);
  newtree->Branch("chi2_x", &br_chi2_x);
  newtree->Branch("chi2_y", &br_chi2_y);
  newtree->Branch("angle_x", &br_angle_x);
  newtree->Branch("angle_y", &br_angle_y);
  newtree->Branch("n_tracks", &br_n_tracks);
  newtree->Branch("total_clusters", &br_total_clusters);
  newtree->Branch("cluster_col", &br_cluster_col);
  newtree->Branch("cluster_row", &br_cluster_row);
  newtree->Branch("cluster_charge", &br_cluster_charge);
  newtree->Branch("cluster_xpos_tel", &br_cluster_xpos_tel);
  newtree->Branch("cluster_ypos_tel", &br_cluster_ypos_tel);
  newtree->Branch("cluster_xpos_local", &br_cluster_xpos_local);
  newtree->Branch("cluster_ypos_local", &br_cluster_ypos_local);
  newtree->Branch("residuals_x", &br_residuals_x);
  newtree->Branch("residuals_y", &br_residuals_y);
  newtree->Branch("residuals", &br_residuals);
  newtree->Branch("cluster_size", &br_cluster_size);
  newtree->Branch("track_x", &br_track_x);
  newtree->Branch("track_y", &br_track_y);
}

void FileWriterTracking::fillTree(){
  newtree->Fill();
}

void FileWriterTracking::saveTree(){

  newfile->cd();
  newtree->Write();
  if (names != nullptr) {
    names->Write(); }
  newfile->Write();
  newfile->Close();
  delete newfile;
}

void FileWriterTracking::clearVectors(){

  for (uint8_t iRoc = 0; iRoc != n_rocs_; iRoc++) {
    br_cluster_col->at(iRoc).clear();
    br_cluster_row->at(iRoc).clear();
    br_cluster_xpos_tel->at(iRoc).clear();
    br_cluster_ypos_tel->at(iRoc).clear();
    br_cluster_xpos_local->at(iRoc).clear();
    br_cluster_ypos_local->at(iRoc).clear();
    br_cluster_charge->at(iRoc).clear();
    br_cluster_size->at(iRoc).clear();
    br_residuals_x->at(iRoc).clear();
    br_residuals_y->at(iRoc).clear();
    br_residuals->at(iRoc).clear();
    br_track_x->at(iRoc).clear();
    br_track_y->at(iRoc).clear();
  }
}

void FileWriterTracking::resizeVectors() {

  br_residuals_x->resize(n_rocs_);
  br_residuals_y->resize(n_rocs_);
  br_residuals->resize(n_rocs_);
  br_track_x->resize(n_rocs_);
  br_track_y->resize(n_rocs_);

  br_cluster_size->resize(n_rocs_);

  br_cluster_col->resize(n_rocs_);
  br_cluster_row->resize(n_rocs_);
  br_cluster_xpos_tel->resize(n_rocs_);
  br_cluster_ypos_tel->resize(n_rocs_);
  br_cluster_xpos_local->resize(n_rocs_);
  br_cluster_ypos_local->resize(n_rocs_);

  br_cluster_charge->resize(n_rocs_);
  is_aligned->resize(n_rocs_, true);
  br_aligned->resize(n_rocs_, true);
}

void FileWriterTracking::setSResidual(uint8_t iRoc, bool def) {
  /** Fill the single cluster residuals. */
  br_sres_x[iRoc] = def ? br_residuals_x->at(iRoc).at(0) : DEF_VAL;
  br_sres_y[iRoc] = def ? br_residuals_y->at(iRoc).at(0) : DEF_VAL;
  br_sres[iRoc] = def ? br_residuals->at(iRoc).at(0) : DEF_VAL;
}

void FileWriterTracking::set_dut_tracks(const vector<float> * z_dut) {
  /** Fill the intersection of the tracks with the DUT planes */
  if (FR_->NTracks() == 1) {
    auto * track = FR_->Track(0);
    for (uint8_t i(0); i < z_dut->size(); i++) {
      float x_pos(track->ExtrapolateX(z_dut->at(i))), y_pos(track->ExtrapolateY(z_dut->at(i)));
      br_dia_track_pos_x[i] = x_pos; br_dia_track_pos_y[i] = y_pos;
      auto pos_loc = FR_->GetAlignment()->TtoLXY(x_pos, y_pos, FR_->Channel(), int(n_rocs_ - n_duts_ + i));
      br_dia_track_pos_x_loc[i] = pos_loc.first; br_dia_track_pos_y_loc[i] = pos_loc.second;
      br_dist_to_dia[i] = sqrt(x_pos * x_pos + y_pos * y_pos);
    }
  } else {
    for (uint8_t i(0); i < n_duts_; i++) {
      br_dia_track_pos_x[i]=br_dia_track_pos_y[i]=br_dia_track_pos_x_loc[i]=br_dia_track_pos_y_loc[i]=DEF_VAL;
    }
  }
}

