#include "FileWriterTracking.h"
#include "TInterpreter.h"

using std::cout; using std::string; using std::stringstream; using std::vector; using std::endl;

/** ============================
 CONSTRUCTOR
 =================================*/
FileWriterTracking::FileWriterTracking(string InFileName, uint8_t telescopeID, PSIFileReader * FR):
  nRoc(GetNumberOfROCS(telescopeID)), nHits(4), TelescopeID(telescopeID) {

  gInterpreter->GenerateDictionary("vector<vector<UShort_t> >");  // add vector<vector> to root dicts
  gInterpreter->GenerateDictionary("vector<vector<Float_t> >");  // add vector<vector> to root dicts
  NewFileName = getFileName(InFileName);
  intree = ((PSIRootFileReader*) FR)->fTree;
  names = ((PSIRootFileReader*) FR)->fMacro;
  newfile = new TFile(NewFileName.c_str(), "RECREATE");
  newtree = intree->CloneTree(0);

  cout << "<<<0>>>>" << endl;
  /** init vectors */
  br_dia_track_pos_x = new vector<float>;
  br_dia_track_pos_y = new vector<float>;
  br_dist_to_dia = new vector<float>;
  br_residuals_x = new vector<vector<float> >;
  br_residuals_y = new vector<vector<float> >;
  br_residuals = new vector<vector<float> >;
  br_single_cluster_residuals = new vector<float>;
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
  cout << "<<<1>>>>" << endl;
  resizeVectors();
  cout << "<<<2>>>>" << endl;
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
  newtree->Branch("total_hits", &br_total_hits);
  newtree->Branch("total_clusters", &br_total_clusters);
  newtree->Branch("n_clusters", &br_n_clusters);
//  newtree->Branch("cluster_plane", &br_cluster_plane);
  newtree->Branch("cluster_col", &br_cluster_col);
  newtree->Branch("cluster_row", &br_cluster_row);
  newtree->Branch("cluster_charge", &br_cluster_charge);
  newtree->Branch("cluster_xpos_tel", &br_cluster_xpos_tel);
  newtree->Branch("cluster_ypos_tel", &br_cluster_ypos_tel);
  newtree->Branch("cluster_xpos_local", &br_cluster_xpos_local);
  newtree->Branch("cluster_ypos_local", &br_cluster_ypos_local);
  newtree->Branch("n_hits", &br_n_hits);
  newtree->Branch("residuals_x", &br_residuals_x);
  newtree->Branch("residuals_y", &br_residuals_y);
  newtree->Branch("residuals", &br_residuals);
  newtree->Branch("s_residuals", &br_single_cluster_residuals);
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
  if (names != nullptr)
    names->Write();
  newfile->Write();
  newfile->Close();
  delete newfile;
}

void FileWriterTracking::clearVectors(){

//  br_cluster_plane.clear();
  for (uint8_t iRoc = 0; iRoc != nRoc; iRoc++) {
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
  br_dia_track_pos_x->clear();
  br_dia_track_pos_y->clear();
  br_dist_to_dia->clear();
}

void FileWriterTracking::resizeVectors() {

  br_residuals_x->resize(nRoc);
  br_residuals_y->resize(nRoc);
  br_residuals->resize(nRoc);
  br_single_cluster_residuals->resize(nRoc);
  br_track_x->resize(nRoc);
  br_track_y->resize(nRoc);

  br_n_clusters.resize(nRoc);
  br_n_hits.resize(nRoc);
  br_cluster_size->resize(nRoc);

  br_cluster_col->resize(nRoc);
  br_cluster_row->resize(nRoc);
  br_cluster_xpos_tel->resize(nRoc);
  br_cluster_ypos_tel->resize(nRoc);
  br_cluster_xpos_local->resize(nRoc);
  br_cluster_ypos_local->resize(nRoc);

  br_cluster_charge->resize(nRoc);
}

