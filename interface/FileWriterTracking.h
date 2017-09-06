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
    uint8_t const nRoc;
    size_t const nHits;

    /** keep track */
    TTree * intree;
    TFile * newfile;
    TTree * newtree;
    TMacro * names;
    std::string NewFileName;

    int16_t TelescopeID;

    /** branch variables */
    uint8_t br_hit_plane_bits;
    std::vector<float> * br_dia_track_pos_x;
    std::vector<float> * br_dia_track_pos_y;
    std::vector<float> * br_dist_to_dia;
    float   br_chi2;
    float   br_chi2_x, br_chi2_y;
    float   br_angle_x, br_angle_y;
    uint8_t br_n_tracks, br_n_clusters;
    std::vector<uint8_t> br_clusters_per_plane;
    std::vector<uint16_t> br_n_hits;
    std::vector<uint16_t> br_cluster_plane;
    std::vector<uint16_t> br_cluster_col;
    std::vector<uint16_t> br_cluster_row;
    std::vector<float> br_cluster_xpos_tel;
    std::vector<float> br_cluster_ypos_tel;
    std::vector<float> br_cluster_xpos_local;
    std::vector<float> br_cluster_ypos_local;
    std::vector<float> br_cluster_charge;
    std::vector<std::vector<float>* > br_charge_all;
    /** double vectors */
    std::vector<std::vector<uint16_t> > * br_cluster_size;
    std::vector<std::vector<float> > * br_residuals_x;
    std::vector<std::vector<float> > * br_residuals_y;
    std::vector<std::vector<float> > * br_residuals;
    float br_coincidence_map;
    std::vector<std::vector<float> *> br_track_x;
    std::vector<std::vector<float> *> br_track_y;
    std::vector<std::vector<float> *> br_smallest_hit_charge;
    std::vector<std::vector<int> *> br_smallest_hit_adc;
    std::vector<std::vector<int> *> br_smallest_hit_pos_col;
    std::vector<std::vector<int> *> br_smallest_hit_pos_row;

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
    uint8_t HitPlaneBits() { return br_hit_plane_bits; }
    uint8_t nTracks() { return br_n_tracks; }
    uint8_t nClusters() { return br_n_clusters; }
    float DiaTrackX(uint8_t roc) { return br_dia_track_pos_x->at(roc); }
    float DiaTrackY(uint8_t roc) { return br_dia_track_pos_y->at(roc); }
    float   AngleX() { return br_angle_x; }
    float   AngleY() { return br_angle_y; }
    float   Chi2() { return br_chi2; }
    float   Chi2X() { return br_chi2_x; }
    float   Chi2Y() { return br_chi2_y; }
    std::vector<std::vector<float>* > ChargeAll() { return br_charge_all; }
    std::vector<std::vector<uint16_t> > * ClusterSize() { return br_cluster_size; }
    size_t GetNHits() { return nHits; }
    float   GetCoincidenceMap() { return br_coincidence_map; }
    std::vector<uint16_t> ClusterRow() { return br_cluster_row; }
    std::vector<uint16_t> ClusterCol() { return br_cluster_col; }
    std::vector<std::vector<float>* > TrackX() { return br_track_x; }
    std::vector<std::vector<float>* > TrackY() { return br_track_y; }
    std::vector<std::vector<float>* > SmallestHitCharge() { return br_smallest_hit_charge; }
    std::vector<std::vector<int>* > SmallestHitADC() { return br_smallest_hit_adc; }
    std::vector<std::vector<int>* > SmallestHitPosCol() { return br_smallest_hit_pos_col; }
    std::vector<std::vector<int>* > SmallestHitPosRow() { return br_smallest_hit_pos_row; }
    std::vector<std::vector<float> > * ResidualX() { return br_residuals_x; }
    std::vector<std::vector<float> > * ResidualY() { return br_residuals_y; }

    /** ============================
        SETTER METHODS
        =================================*/
    void setHitPlaneBits(uint8_t value) { br_hit_plane_bits = value; }
    void setNTracks(uint8_t value) { br_n_tracks = value; }
    void setNClusters(uint8_t value) { br_n_clusters = value; }
    void setDiaTracks(float xVal, float yVal) { br_dia_track_pos_x->push_back(xVal); br_dia_track_pos_y->push_back(yVal); }
    void setDistDia(float xVal, float yVal) { br_dist_to_dia->push_back(sqrt(xVal * xVal + yVal * yVal)); }
    void setAngleX(float value) { br_angle_x = value; }
    void setAngleY(float value) { br_angle_y = value; }
    void setChi2(float value) { br_chi2 = value; }
    void setChi2X(float value) { br_chi2_x = value; }
    void setChi2Y(float value) { br_chi2_y = value; }
    void setChargeAll(uint8_t iRoc, float value) { br_charge_all[iRoc]->push_back(value); }
    void setClusterSize(uint8_t iRoc, int value) { br_cluster_size->at(iRoc).push_back(value); }
    void setClusters(uint8_t iRoc, uint8_t value) { br_clusters_per_plane[iRoc] = value; }
    void setClusterPlane(uint16_t value) { br_cluster_plane.push_back(value); }
    void setNHits(uint8_t iRoc, uint16_t value) { br_n_hits[iRoc] = value; }
    void setCoincidenceMap(float value) { br_coincidence_map = value; }
    void setClusterColumn(int value) { br_cluster_col.push_back(value); }
    void setClusterRow(int value) { br_cluster_row.push_back(value); }
    void setClusterXPosTel(float value) { br_cluster_xpos_tel.push_back(value); }
    void setClusterYPosTel(float value) { br_cluster_ypos_tel.push_back(value); }
    void setClusterXPosLocal(float value) { br_cluster_xpos_local.push_back(value); }
    void setClusterYPosLocal(float value) { br_cluster_ypos_local.push_back(value); }
    void setClusterCharge(float value) { br_cluster_charge.push_back(value); }
    void setTrackX(uint8_t iRoc, float value) { br_track_x[iRoc]->push_back(value); }
    void setTrackY(uint8_t iRoc, float value) { br_track_y[iRoc]->push_back(value); }
    void setSmallestHitCharge(uint8_t iRoc, float value) { br_smallest_hit_charge[iRoc]->push_back(value); }
    void setSmallestHitCharge(uint8_t iRoc, int value) { br_smallest_hit_adc[iRoc]->push_back(value); }
    void setSmallestHitPosCol(uint8_t iRoc, int value) { br_smallest_hit_pos_col[iRoc]->push_back(value); }
    void setSmallestHitPosRow(uint8_t iRoc, int value) { br_smallest_hit_pos_row[iRoc]->push_back(value); }
    void setResidualX(uint8_t iRoc, float value) { br_residuals_x->at(iRoc).push_back(value); }
    void setResidualY(uint8_t iRoc, float value) { br_residuals_y->at(iRoc).push_back(value); }
    void setResidual(uint8_t iRoc, float value) { br_residuals->at(iRoc).push_back(value); }

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
