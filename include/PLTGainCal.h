#ifndef GUARD_PLTGainCal_h
#define GUARD_PLTGainCal_h

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>

#include "TString.h"
#include "TMath.h"
#include "TF1.h"

#include "PLTHit.h"
#include "PLTCluster.h"
#include "PLTPlane.h"
#include "PLTU.h"
#include "GetNames.h"


class PLTGainCal
{
  public:
    PLTGainCal ();
    PLTGainCal (int, bool); // number of ROCs, isExternalFunction
    PLTGainCal (std::string const, int const NParams=5);
    ~PLTGainCal ();


    static int const DEBUGLEVEL = 0;

    void SetCharge ( PLTHit&, uint8_t telescopeID=10);
    float GetCharge(int const ch, int const telescopeID, int const roc, int const col, int const row, int adc);

    void ReadGainCalFile (std::string const GainCalFileName, int roc);
    void ReadGainCalFile3 (std::string const GainCalFileName);
    void ReadGainCalFile5 (std::string const GainCalFileName);
    void ReadGainCalFileExt (std::string const GainCalFileName, int const roc=0);
    void ReadTesterGainCalFile (std::string const GainCalFileName);
    void ReadVcalCal();

    void CheckGainCalFile (std::string const GainCalFileName, int const Channel);

    void PrintGainCal5 ();
    void PrintGainCal (FILE* f = 0x0);

    static int RowIndex (int const);
    static int ColIndex (int const);
    static int ChIndex (int const);
    static int RocIndex (int const);

    float GetCoef(int const, int const, int const, int const, int const);

    bool IsGood () { return fIsGood; }

    int GetHardwareID (int const);

    void ResetGC ();



  private:
    bool fIsGood;
    bool fIsExternalFunction;

    int  fNParams {}; // how many parameters for this gaincal
    TF1 fFitFunction;

    static int const MAXCHNS =   1;
    static int const MAXROWS =  80;
    static int const MAXCOLS =  52;
    static int const MAXROCS =   6;

    static int const NCHNS =   1;

    int const NROCS;

    // Switched from
    // float GC[NCHNS][NROCS][NCOLS][NROWS][6]
    std::vector<std::vector<std::vector<std::vector<std::vector<float > > > > > GC;
    std::vector<std::pair<float, float>> VC;  // VCal Calibration


    // Map for hardware locations by fed channel
    std::map<int, int> fHardwareMap;

};















#endif
