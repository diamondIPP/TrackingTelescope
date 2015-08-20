#include "GetNames.h"

using namespace std;

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
        else if (telescopeID==10)    return "ALIGNMENT/Alignment_ETHTelescope_telescope10.dat";
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
    else if (telescopeID == 10)  return "outer_pixel_masks/outerPixelMask_Telescope10.txt";
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
    else if (telescopeID == -1) return "calibration_lists/GKCalibrationList_Telescope5.txt";
    else {
        cout << "ERROR: No Calibration file for telescopeID=" << telescopeID << endl;
        cout << "Exiting.." << endl;
        exit(0);
    }
}

int GetNumberOfROCS(int telescopeID){

    if ((telescopeID == 1) || (telescopeID == 2) || (telescopeID == 3) || (telescopeID == 8))
        return 6;
    else if (telescopeID == 4)
        return 2;
    else if (telescopeID == 10)
        return 7;
    else if ((telescopeID == 5) || (telescopeID == 6) || (telescopeID == 7) || (telescopeID == -1) || (telescopeID == 9))
        return 4;
    else {
        cout << "ERROR: Number of ROCs not defined for telescopeID=" << telescopeID << endl;
        cout << "Exiting.." << endl;
        exit(0);
    }
}

int GetUseGainInterpolator(int telescopeID){

    if (telescopeID == 2)   return true;
    else                    return false;
}

int GetUseExternalCalibrationFunction(int telescopeID){

    if (telescopeID == 7 || telescopeID == 8 || telescopeID == 9 || telescopeID == 10)   return true;
    else                    return false;
}

int GetUseRootInput(int telescopeID){

  if ((telescopeID == -1) || (telescopeID == 7) || telescopeID == 8 || telescopeID == 9 || telescopeID == 10)
        return true;
  else  return false;
}

int GetNumberOfSignals(int telescopeID){

    uint8_t id = telescopeID;
    if ((id == 7) || (id == 8) || (id == 9) || (id == 10))
        return 4;
    else {
        cerr << "ERROR: Number of Signals is not defined for telescopeID=" << telescopeID << endl;
        cout << "Exiting.." << endl;
        exit(0);
    }
}

const char * GetSignalBranchName(){

    return "sig";
}

