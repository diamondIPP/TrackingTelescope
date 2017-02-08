#ifndef FileWriterTracking_h
#define FileWriterTracking_h

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <string>

#include "TTree.h"
#include "TFile.h"

#include "PSIRootFileReader.h"

using namespace std;

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
    string NewFileName;

    int16_t TelescopeID;

    /** branch variables */
    uint8_t br_hit_plane_bits;
    float   br_diam1_track_x, br_diam1_track_y;
    float   br_diam2_track_x, br_diam2_track_y;
    float   br_dist_to_dia1, br_dist_to_dia2;
    float   br_chi2;
    float   br_chi2_x, br_chi2_y;
    float   br_angle_x, br_angle_y;
    float   br_slope_x, br_slope_y;
    uint8_t br_n_tracks, br_n_clusters;
    vector<uint8_t> br_clusters_per_plane;
    vector<uint16_t> br_n_hits;
    vector<uint16_t> br_cluster_plane;
    vector<uint16_t> br_cluster_col;
    vector<uint16_t> br_cluster_row;
    vector<float> br_cluster_xpos_tel;
    vector<float> br_cluster_ypos_tel;
    vector<float> br_cluster_xpos_local;
    vector<float> br_cluster_ypos_local;
    vector<float> br_cluster_charge;
    vector<vector<float>* > br_charge_all;
    vector<vector<int>* > br_cluster_size;
    float br_coincidence_map;
    vector<vector<float> *> br_track_x;
    vector<vector<float> *> br_track_y;
    vector<vector<float> *> br_smallest_hit_charge;
    vector<vector<int> *> br_smallest_hit_adc;
    vector<vector<int> *> br_smallest_hit_pos_col;
    vector<vector<int> *> br_smallest_hit_pos_row;
    vector<vector<float> *> br_residual_local_x;
    vector<vector<float> *> br_residual_local_y;
    vector<float> br_residuals_x;
    vector<float> br_residuals_y;

    /** some functions*/
    string getFileName(string);

public:

    /** ============================
     CONSTRUCTOR
     =================================*/
    FileWriterTracking(string, uint8_t, PSIFileReader * FR);
    ~FileWriterTracking();


    /** ============================
     GET-FUNCTIONS
     =================================*/
    string FileName() {return NewFileName; }
    TTree * InTree() { return intree; }
    uint8_t HitPlaneBits() { return br_hit_plane_bits; }
    uint8_t nTracks() { return br_n_tracks; }
    uint8_t nClusters() { return br_n_clusters; }
    float   Dia1TrackX() { return br_diam1_track_x; }
    float   Dia1TrackY() { return br_diam1_track_y; }
    float   Dia2TrackX() { return br_diam2_track_x; }
    float   Dia2TrackY() { return br_diam2_track_y; }
    float   AngleX() { return br_angle_x; }
    float   AngleY() { return br_angle_y; }
    float   SlopeX() { return br_slope_x; }
    float   SlopeY() { return br_slope_y; }
    float   Chi2() { return br_chi2; }
    float   Chi2X() { return br_chi2_x; }
    float   Chi2Y() { return br_chi2_y; }
    vector<vector<float>* > ChargeAll() { return br_charge_all; }
    vector<vector<int>* > ClusterSize() { return br_cluster_size; }
    size_t GetNHits() { return nHits; }
    float   GetCoincidenceMap() { return br_coincidence_map; }
    vector<uint16_t> ClusterRow() { return br_cluster_row; }
    vector<uint16_t> ClusterCol() { return br_cluster_col; }
    vector<vector<float>* > TrackX() { return br_track_x; }
    vector<vector<float>* > TrackY() { return br_track_y; }
    vector<vector<float>* > SmallestHitCharge() { return br_smallest_hit_charge; }
    vector<vector<int>* > SmallestHitADC() { return br_smallest_hit_adc; }
    vector<vector<int>* > SmallestHitPosCol() { return br_smallest_hit_pos_col; }
    vector<vector<int>* > SmallestHitPosRow() { return br_smallest_hit_pos_row; }
    vector<vector<float>* > ResidualLocalX() { return br_residual_local_x; }
    vector<vector<float>* > ResidualLocalY() { return br_residual_local_y; }
    vector<float> ResidualX() { return br_residuals_x; }
    vector<float> ResidualY() { return br_residuals_y; }

    /** ============================
     SET-FUNCTIONS
     =================================*/
    void setHitPlaneBits(uint8_t value) { br_hit_plane_bits = value; }
    void setNTracks(uint8_t value) { br_n_tracks = value; }
    void setNClusters(uint8_t value) { br_n_clusters = value; }
    void setDia1TrackX(float value) { br_diam1_track_x = value; }
    void setDia1TrackY(float value) { br_diam1_track_y = value; }
    void setDia2TrackX(float value) { br_diam2_track_x = value; }
    void setDia2TrackY(float value) { br_diam2_track_y = value; }
    void setDistDia1(float xVal, float yVal) { br_dist_to_dia1 = sqrt(xVal*xVal + yVal*yVal); }
    void setDistDia2(float xVal, float yVal) { br_dist_to_dia2 = sqrt(xVal*xVal + yVal*yVal); }
    void setAngleX(float value) { br_angle_x = value; }
    void setAngleY(float value) { br_angle_y = value; }
    void setSlopeX(float value) { br_slope_x = value; }
    void setSlopeY(float value) { br_slope_y = value; }
    void setChi2(float value) { br_chi2 = value; }
    void setChi2X(float value) { br_chi2_x = value; }
    void setChi2Y(float value) { br_chi2_y = value; }
    void setChargeAll(uint8_t iRoc, float value) { br_charge_all[iRoc]->push_back(value); }
    void setClusterSize(uint8_t iRoc, int value) { br_cluster_size[iRoc]->push_back(value); }
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
    void setResidualLocalX(uint8_t iRoc, float value) { br_residual_local_x[iRoc]->push_back(value); }
    void setResidualLocalY(uint8_t iRoc, float value) { br_residual_local_y[iRoc]->push_back(value); }
    void setResidualsX(uint8_t iRoc, float value) { br_residuals_x.at(iRoc) = value; }
    void setResidualsY(uint8_t iRoc, float value) { br_residuals_y.at(iRoc) = value; }

    /** ============================
     AUXILIARY FUNCTIONS
     =================================*/
    void addBranches();
    void saveTree();
    void fillTree();
    void clearVectors();


};



#endif // FileWriterTracking_h
