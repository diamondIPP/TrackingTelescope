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
    float   br_chi2;
    float   br_chi2_x, br_chi2_y;
    float   br_slope_x, br_slope_y;
    uint8_t br_n_tracks, br_n_clusters;
    vector<vector<float>* > br_charge_all;
//    br_charge_all.resize(NROC);

    /** some functions*/
    string getFileName(string);

public:

    /** ============================
     CONSTRUCTOR
     =================================*/
    FileWriterTracking(string, uint8_t, PSIFileReader * FR);
    ~FileWriterTracking();


    string FileName() {return NewFileName; }
    TTree * InTree() { return intree; }





};



#endif // FileWriterTracking_h
