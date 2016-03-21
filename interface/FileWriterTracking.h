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

    /** keep track */
    TTree * intree;
    TFile * newfile;
    TTree * newtree;
    string NewFileName;

    /** branch variables */
    uint8_t br_hit_plane_bits;
    float   br_diam1_track_x, br_diam1_track_y;
    float   br_diam2_track_x, br_diam2_track_y;
    float   br_dist_to_dia1, br_dist_to_dia2;
    float   br_chi2;
    float   br_chi2_x, br_chi2_y;
    float   br_slope_x, br_slope_y;
    uint8_t br_n_tracks, br_n_clusters;
    vector<uint8_t> br_clusters_per_plane;
    vector<vector<float>* > br_charge_all;
    // Add vector for pulse height of each event corresponding to each adc - DA
    std::vector<int> *br_pulse_height;
//    vector<vector<float> > br_cluster_pos_x;
//    vector<vector<float> > br_cluster_pos_y;
//    vector<vector<float> > br_test;

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
    float   SlopeX() { return br_slope_x; }
    float   SlopeY() { return br_slope_y; }
    float   Chi2() { return br_chi2; }
    float   Chi2X() { return br_chi2_x; }
    float   Chi2Y() { return br_chi2_y; }
    vector<vector<float>* > ChargeAll() { return br_charge_all; }

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
    void setSlopeX(float value) { br_slope_x = value; }
    void setSlopeY(float value) { br_slope_y = value; }
    void setChi2(float value) { br_chi2 = value; }
    void setChi2X(float value) { br_chi2_x = value; }
    void setChi2Y(float value) { br_chi2_y = value; }
    void setChargeAll(uint8_t iRoc, float value) { br_charge_all[iRoc]->push_back(value); }
    void setClusters(uint8_t iRoc, uint8_t value) { br_clusters_per_plane[iRoc] = value; }
//    void setClusterPositionX(uint8_t iRoc, float value) { br_cluster_pos_x[iRoc].push_back(value); }
//    void setClusterPositionY(uint8_t iRoc, float value) { br_cluster_pos_y[iRoc].push_back(value); }


    /** ============================
     AUXILIARY FUNCTIONS
     =================================*/
    void addBranches();
    void saveTree();
    void fillTree();
    void clearVectors();


};



#endif // FileWriterTracking_h
