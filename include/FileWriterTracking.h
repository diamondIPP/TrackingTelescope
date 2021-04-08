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

class FileWriterTracking{

private:
    /** constants */
    const uint16_t nRoc;
    const size_t nHits;

    /** keep track */
    TTree * intree;
    TFile * newfile;
    TTree * newtree;
    TMacro * names;
    std::string NewFileName;

    int16_t TelescopeID;

    /** ============================
        BRANCH VARIABLES
     =================================*/
    uint16_t br_hit_plane_bits;
    std::vector<bool> * is_aligned;
    std::vector<bool> * br_aligned;

    /** tracks */
    uint8_t br_n_tracks;
    std::vector<float> * br_dia_track_pos_x, * br_dia_track_pos_y;
    std::vector<float> * br_dia_track_pos_x_loc, * br_dia_track_pos_y_loc;
    std::vector<float> * br_dist_to_dia;
    float   br_chi2, br_chi2_x, br_chi2_y;
    float   br_angle_x, br_angle_y;
    std::vector<std::vector<float> > * br_residuals_x, * br_residuals_y;
    std::vector<std::vector<float> > * br_residuals;
    std::vector<float> * br_single_cluster_residuals;
    std::vector<std::vector<float> > * br_track_x, * br_track_y;

    /** cluster numbers */
    std::vector<uint16_t> br_n_hits;
    uint8_t br_total_clusters;
    std::vector<uint8_t> br_n_clusters;
    std::vector<std::vector<uint16_t> > * br_cluster_size;
//    std::vector<uint16_t> br_cluster_plane;

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
    FileWriterTracking(std::string, uint8_t, PSIFileReader * FR);
    ~FileWriterTracking();


    /** ============================
     GET-FUNCTIONS
     =================================*/
    std::string FileName() {return NewFileName; }
    TTree * InTree() { return intree; }
    uint16_t HitPlaneBits() { return br_hit_plane_bits; }
    uint8_t nTracks() { return br_n_tracks; }
    uint8_t nClusters() { return br_total_clusters; }
    float DiaTrackX(uint8_t roc) { return br_dia_track_pos_x->at(roc); }
    float DiaTrackY(uint8_t roc) { return br_dia_track_pos_y->at(roc); }
    float   AngleX() { return br_angle_x; }
    float   AngleY() { return br_angle_y; }
    float   Chi2() { return br_chi2; }
    float   Chi2X() { return br_chi2_x; }
    float   Chi2Y() { return br_chi2_y; }
    std::vector<std::vector<uint16_t> > * ClusterSize() { return br_cluster_size; }
    std::vector<uint16_t> GetNHits() { return br_n_hits; }
    std::vector<std::vector<float> > * TrackX() { return br_track_x; }
    std::vector<std::vector<float> > * TrackY() { return br_track_y; }
    std::vector<std::vector<float> > * ResidualX() { return br_residuals_x; }
    std::vector<std::vector<float> > * ResidualY() { return br_residuals_y; }
    bool lastIsAligned(uint8_t iRoc) { return is_aligned->at(iRoc); }

    /** ============================
        SETTER METHODS
        =================================*/
    void setHitPlaneBits(uint16_t value) { br_hit_plane_bits = value; }
    /** tracks */
    void setNTracks(uint8_t value) { br_n_tracks = value; }
    void setDiaTracks(float xVal, float yVal) { br_dia_track_pos_x->push_back(xVal); br_dia_track_pos_y->push_back(yVal); }
    void setDiaTracksLocal(std::pair<float, float> pos) { br_dia_track_pos_x_loc->push_back(pos.first); br_dia_track_pos_y_loc->push_back(pos.second); }
    void setDistDia(float xVal, float yVal) { br_dist_to_dia->push_back(sqrt(xVal * xVal + yVal * yVal)); }
    void setAngle(float xval, float yval) { br_angle_x = xval;  br_angle_y = yval; }
    void setChi2(float total, float xval, float yval) { br_chi2 = total; br_chi2_x = xval; br_chi2_y = yval; }
    void setResidualXY(uint8_t iRoc, float x, float y) { br_residuals_x->at(iRoc).push_back(x); br_residuals_y->at(iRoc).push_back(y); }
    void setResidual(uint8_t iRoc, float value) { br_residuals->at(iRoc).push_back(value); }
    void setSResidual(uint8_t iRoc, float value) { br_single_cluster_residuals->at(iRoc) = value; }
    void setTrackPos(uint8_t iRoc, float x, float y) { br_track_x->at(iRoc).push_back(x); br_track_x->at(iRoc).push_back(y); }
    /** cluster numbers */
    void setTotalClusters(uint8_t value) { br_total_clusters = value; }
    void setNHits(uint8_t iRoc, uint16_t value) { br_n_hits.at(iRoc) = value; }
    void setNClusters(uint8_t iRoc, uint8_t value) { br_n_clusters.at(iRoc) = value; }
    void setClusterSize(uint8_t iRoc, int value) { br_cluster_size->at(iRoc).push_back(value); }
    /** cluster positions */
    void setClusterPos(uint8_t iRoc, int col, int row) { br_cluster_col->at(iRoc).push_back(col); br_cluster_row->at(iRoc).push_back(row);}
    void setClusterPosLocal(uint8_t iRoc, float x, float y) { br_cluster_xpos_local->at(iRoc).push_back(x); br_cluster_ypos_local->at(iRoc).push_back(y); }
    void setClusterPosTel(uint8_t iRoc, float x, float y) { br_cluster_xpos_tel->at(iRoc).push_back(x); br_cluster_ypos_tel->at(iRoc).push_back(y); }
    /** cluster charge */
    void setClusterCharge(uint8_t iRoc, float value) { br_cluster_charge->at(iRoc).push_back(value); }
//    void setClusterPlane(uint16_t value) { br_cluster_plane.push_back(value); }
    void setAligned(uint8_t iRoc, bool aligned) { br_aligned->at(iRoc) = aligned; }
    void setOldAligned(uint8_t iRoc, bool aligned) { is_aligned->at(iRoc) = aligned; }

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
