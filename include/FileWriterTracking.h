#ifndef FileWriterTracking_h
#define FileWriterTracking_h

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>
#include <map>

#include "TTree.h"
#include "TFile.h"

#include "PSIRootFileReader.h"

#define DEF_VAL -999

class FileWriterTracking{

private:
    PSIFileReader * FR_;

    /** constants */
    const uint16_t n_rocs_;
    const uint8_t n_duts_;

    /** keep track */
    TTree * intree;
    TFile * newfile;
    TTree * newtree;
    TMacro * names;
    std::string NewFileName;

    /** ============================
        BRANCH VARIABLES
     =================================*/
    uint16_t br_hit_plane_bits = 0;
    std::vector<bool> * is_aligned;
    std::vector<bool> * br_aligned;

    /** tracks */
    uint8_t br_n_tracks = 0;
    float *br_dia_track_pos_x, *br_dia_track_pos_y, *br_dia_track_pos_x_loc, *br_dia_track_pos_y_loc, *br_dist_to_dia;
    float br_chi2 = DEF_VAL, br_chi2_x = DEF_VAL, br_chi2_y = DEF_VAL;
    float br_angle_x = DEF_VAL, br_angle_y = DEF_VAL;
    std::vector<std::vector<float> > * br_residuals_x, * br_residuals_y;
    std::vector<std::vector<float> > * br_residuals;
    float *br_sres_x, *br_sres_y, *br_sres;
    std::vector<std::vector<float> > * br_track_x, * br_track_y;

  /** cluster numbers */
    uint8_t br_total_clusters = 0;
    uint16_t * br_n_hits;
    uint8_t * br_n_clusters;
    std::vector<std::vector<uint16_t> > * br_cluster_size;

    /** cluster positions */
    std::vector<std::vector<uint16_t> > * br_cluster_col, * br_cluster_row;
    std::vector<std::vector<float> > * br_cluster_xpos_tel, * br_cluster_ypos_tel;      // telescope coordinates (after alignment)
    std::vector<std::vector<float> > * br_cluster_xpos_local, * br_cluster_ypos_local;  // local coordinates

    /** cluster charge */
    std::vector<std::vector<float> > * br_cluster_charge;

    /** some functions*/
    std::string getFileName(std::string &);

public:

    /** ============================
     CONSTRUCTOR
     =================================*/
    FileWriterTracking(std::string, PSIFileReader * FR);


    /** ============================
     GET-FUNCTIONS
     =================================*/
    TTree * InTree() { return intree; }
    uint8_t nClusters() const { return br_total_clusters; }

    /** ============================
        SETTER METHODS
        =================================*/
    void setHitPlaneBits(uint16_t value) { br_hit_plane_bits = value; }
    /** tracks */
    void setNTracks(uint8_t value) { br_n_tracks = value; }
    void set_dut_tracks(const std::vector<float>*);
    void setAngle(float xval, float yval) { br_angle_x = xval;  br_angle_y = yval; }
    void setChi2(float total, float xval, float yval) { br_chi2 = total; br_chi2_x = xval; br_chi2_y = yval; }
    void setResidualXY(uint8_t iRoc, float x, float y) { br_residuals_x->at(iRoc).push_back(x); br_residuals_y->at(iRoc).push_back(y); }
    void setResidual(uint8_t iRoc, float value) { br_residuals->at(iRoc).push_back(value); }
    void setSResidual(uint8_t iRoc, bool def);
    void setTrackPos(uint8_t iRoc, float x, float y) { br_track_x->at(iRoc).push_back(x); br_track_x->at(iRoc).push_back(y); }
    /** cluster numbers */
    void setTotalClusters(uint8_t value) { br_total_clusters = value; }
    void setNHits(uint8_t i_roc, uint16_t value) { br_n_hits[i_roc] = value; }
    void setNClusters(uint8_t i_roc, uint8_t value) { br_n_clusters[i_roc] = value; }
    void setClusterSize(uint8_t iRoc, int value) { br_cluster_size->at(iRoc).push_back(value); }
    /** cluster positions */
    void setClusterPos(uint8_t iRoc, int col, int row) { br_cluster_col->at(iRoc).push_back(col); br_cluster_row->at(iRoc).push_back(row);}
    void setClusterPosLocal(uint8_t iRoc, float x, float y) { br_cluster_xpos_local->at(iRoc).push_back(x); br_cluster_ypos_local->at(iRoc).push_back(y); }
    void setClusterPosTel(uint8_t iRoc, float x, float y) { br_cluster_xpos_tel->at(iRoc).push_back(x); br_cluster_ypos_tel->at(iRoc).push_back(y); }
    /** cluster charge */
    void setClusterCharge(uint8_t iRoc, float value) { br_cluster_charge->at(iRoc).push_back(value); }

    /** ============================
     AUXILIARY FUNCTIONS
     =================================*/
    void addBranches();
    void resizeVectors();
    void saveTree();
    void fillTree();
    void clearVectors();

};

#endif // FileWriterTracking_h
