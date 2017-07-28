#include <PLTU.h>
#include "GetNames.h"

using namespace std;

uint16_t nTelescopes = 29;

/** Get the correct alignment for a given telescope */
string GetAlignmentFilename(int telescopeID, bool useInitial){

    /** Initial Alignment (start values for finding alignment) */
    if (useInitial){
        if ((telescopeID==1) || (telescopeID==2) ){
            return "ALIGNMENT/Alignment_ETHTelescope_initial.dat";
        }
        else if (telescopeID==7){
            return "ALIGNMENT/Alignment_ETHTelescope_initial_telescope7.dat";
        }
        else if (telescopeID==10){
            return "ALIGNMENT/Alignment_ETHTelescope_initial_telescope10.dat";
        }
        else if ((telescopeID==5) || (telescopeID==6) || (telescopeID==-1)){
            return "ALIGNMENT/Alignment_ETHTelescope_initial_4planes.dat";
        }
        else {
              cout << "ERROR: No Initial-Alignment file for telescopeID=" << telescopeID << endl;
              cout << "Exiting..." << endl;
              exit(0);
        }
    }

    /** Real Alignment */
    else{
        if (telescopeID==1)         return "ALIGNMENT/Alignment_ETHTelescope_run316.dat";
        else if (telescopeID==2)    return "ALIGNMENT/Alignment_ETHTelescope_run466.dat";
        else if (telescopeID==5)    return "ALIGNMENT/Alignment_ETHTelescope_4planes_run63.dat";
        else if (telescopeID==6)    return "ALIGNMENT/Alignment_ETHTelescope_4planesCERN_run71.dat";
        else if (telescopeID==7)    return "ALIGNMENT/Alignment_ETHTelescope_telescope7.dat";
        else if (telescopeID==8)    return "ALIGNMENT/Alignment_ETHTelescope_telescope8.dat";
        else if (telescopeID==9)    return "ALIGNMENT/Alignment_ETHTelescope_telescope9.dat";
        else if (telescopeID >= 10)   return string(TString::Format("ALIGNMENT/telescope%d.dat", telescopeID));
        else if (telescopeID==-1)   return "ALIGNMENT/Alignment_ETHTelescope_initial_4planes.dat";

        else{
            cout << "ERROR: No Alignment file for telescopeID=" << telescopeID << endl;
            cout << "Exiting..." << endl;
            exit(0);
        }
    }
}

string GetMaskingFilename(int telescopeID){

    if      (telescopeID == 1)  return "outer_pixel_masks/outerPixelMask_Telescope1.txt";
    else if (telescopeID == 2)  return "outer_pixel_masks/outerPixelMask_Telescope2.txt";
    else if (telescopeID == 5)  return "outer_pixel_masks/outerPixelMask_Telescope5.txt";
    else if (telescopeID == 6)  return "outer_pixel_masks/outerPixelMask_Telescope6.txt";
    else if (telescopeID == 7)  return "outer_pixel_masks/outerPixelMask_Telescope7.txt";
    else if (telescopeID == 8)  return "outer_pixel_masks/outerPixelMask_Telescope8.txt";
    else if (telescopeID == 9)  return "outer_pixel_masks/outerPixelMask_Telescope9.txt";
    else if (telescopeID == 22)  return "outer_pixel_masks/outerPixelMask_Telescope22.txt";
    else if (telescopeID >= 25)  return "outer_pixel_masks/outerPixelMask_Telescope22.txt";
    else if (telescopeID >= 10)  return "outer_pixel_masks/outerPixelMask_Telescope10.txt";
    else if (telescopeID == -1) return "outer_pixel_masks/outerPixelMask_Telescope5.txt";
    else {
        cout << "ERROR: No Masking file for telescopeID=" << telescopeID << endl;
        cout << "Exiting.." << endl;
        exit(0);
    }
}

string GetCalibrationFilename(int telescopeID){

    if      (telescopeID == 1)  return "calibration_lists/GKCalibrationList.txt";
    else if (telescopeID == 2)  return "calibration_lists/GKCalibrationList_Telescope2.txt";
    else if (telescopeID == 5)  return "calibration_lists/GKCalibrationList_Telescope5.txt";
    else if (telescopeID == 6)  return "calibration_lists/GKCalibrationList_Telescope6.txt";
    else if (telescopeID == 7)  return "calibration_lists/GKCalibrationList_Telescope7.txt";
    else if (telescopeID == 8)  return "calibration_lists/GKCalibrationList_Telescope8.txt";
    else if (telescopeID == 9)  return "calibration_lists/GKCalibrationList_Telescope9.txt";
    else if (telescopeID == 10)  return "calibration_lists/GKCalibrationList_Telescope10.txt";
    else if (telescopeID == 11)  return "calibration_lists/GKCalibrationList_Telescope9.txt";
    else if (telescopeID == 12)  return "calibration_lists/GKCalibrationList_Telescope12.txt";
    else if (telescopeID == 13)  return "calibration_lists/GKCalibrationList_Telescope13.txt";
    else if (telescopeID == 14)  return "calibration_lists/GKCalibrationList_Telescope12.txt";
    else if (telescopeID == 15)  return "calibration_lists/GKCalibrationList_Telescope13.txt";
    else if (telescopeID == 16)  return "calibration_lists/GKCalibrationList_Telescope7.txt";
    else if (telescopeID == 17)  return "calibration_lists/GKCalibrationList_Telescope12.txt";
    else if (telescopeID == 18)  return "calibration_lists/GKCalibrationList_Telescope12.txt";
    else if (telescopeID == 19)  return "calibration_lists/GKCalibrationList_Telescope12.txt";
    else if (telescopeID == 20)  return "calibration_lists/GKCalibrationList_Telescope12.txt";
    else if (telescopeID == 22)  return "calibration_lists/GKCalibrationList_Telescope22.txt";
    else if (telescopeID == 25)  return "calibration_lists/GKCalibrationList_Telescope25.txt";
    else if (telescopeID >= 10)  return "calibration_lists/GKCalibrationList_Telescope12.txt";
    else if (telescopeID == -1) return "calibration_lists/GKCalibrationList_Telescope5.txt";
    else {
        cout << "ERROR: No Calibration file for telescopeID=" << telescopeID << endl;
        cout << "Exiting.." << endl;
        exit(0);
    }
}

uint8_t GetNumberOfROCS(int16_t telescopeID){


    int16_t id = telescopeID;
    if ((id == 1) || (id == 2) || (id == 3) || (id == 8))
        return 6;
    else if (id == 4)
        return 2;
    else if (id == 10 || id == 13 || id == 15)
        return 7;
    else if(id == 22 or id == 25)
        return 6;
    else if ((id == 5) || (id == 6) || (id == 7) || (id == -1) || (id >= 9))
        return 4;
    else {
        cout << "ERROR: Number of ROCs not defined for telescopeID=" << telescopeID << endl;
        cout << "Exiting.." << endl;
        exit(0);
    }
}

int GetUseGainInterpolator(int telescopeID){

    return telescopeID == 2;
}

int GetUseExternalCalibrationFunction(int telescopeID){

    return telescopeID == 7 || telescopeID == 8 || telescopeID == 9 || telescopeID == 10 || telescopeID >= 11;
}

int GetUseRootInput(int telescopeID){

    return (telescopeID == -1) || (telescopeID == 7) || telescopeID == 8 || telescopeID == 9 || telescopeID == 10 || telescopeID >= 11;
}

bool UseFileWriter(uint8_t telescopeID){

    vector<uint8_t> ids;
    for (uint8_t id = 7; id <= nTelescopes; id++)
        ids.push_back(id);
    return in(telescopeID, ids);
}

bool FillSignalHistos(uint8_t telescopeID){

    vector<uint8_t> ids = {7, 8, 9};
    return in(telescopeID, ids);
}

bool UseDigitalCalibration(int16_t telescopeID){

    vector<int16_t> ids = {10, 13, 15, 22, 25};
    return in(telescopeID, ids);

}

int GetNumberOfSignals(int16_t telescopeID){

    int16_t id = telescopeID;
    if (id == 10 || id == 13 || id == 15 || id == 22 or id == 25)
        return 0;
    else if ((id == 7) || (id == 8) || (id == 9) || id >= 11)
        return 4;
    else {
        cerr << "ERROR: Number of Signals is not defined for telescopeID=" << telescopeID << endl;
        cout << "Exiting.." << endl;
        exit(0);
    }
}

const char * GetSignalBranchName(){

    return "blub";
}

vector<string> & split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

template<typename T>
bool in(T num, vector<T> ids){

    return find(ids.begin(), ids.end(), num ) != ids.end();
}

float GetDiamondZPosition(int16_t id, uint8_t diamond){

    uint8_t tel = 0;
    if (id > 16 && id < 22)
        tel = 1;
    else if (id > 22)
        tel = 2;
    if (diamond == 1)
        return PLTU::DIA1Z[tel];
    else if (diamond == 2)
        return PLTU::DIA2Z[tel];
    else
        return -1;
}
