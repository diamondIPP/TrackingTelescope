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


class PLTGainCal {

public:
  PLTGainCal ();
  PLTGainCal (int, bool); // number of ROCs, isExternalFunction
  PLTGainCal (const std::string&, int);
  ~PLTGainCal () = default;

  static int const DEBUGLEVEL = 0;

  void SetCharge(PLTHit &Hit) { Hit.SetCharge(GetCharge(Hit.Channel(), Hit.ROC(), Hit.Column(), Hit.Row(), Hit.ADC())); }
  float GetCharge(int ch, int roc, int col, int row, int adc);

  void ReadGainCalFile (const std::string & GainCalFileName, int=0);
  void ReadGainCalFile3 (const std::string & GainCalFileName);
  void ReadGainCalFile5 (const std::string & GainCalFileName);
  void ReadGainCalFileExt (const std::string & GainCalFileName, int roc=0);
  void ReadTesterGainCalFile (std::string const GainCalFileName);
  void ReadVcalCal();

  void CheckGainCalFile (std::string const GainCalFileName, int const Channel);
  void PrintGainCal5 ();
  void PrintGainCal (FILE* f = 0x0);

  static int RowIndex (const int row) { return row - PLTU::FIRSTROW; }
  static int ColIndex (const int col) { return col - PLTU::FIRSTCOL; }
  static int ChIndex (const int ch) { return ch - 1; }
  static int RocIndex (const int roc) { return roc; }

  bool IsGood () { return fIsGood; }
  int GetHardwareID (int const);

  void ResetGC ();


private:
  bool fIsGood = false;
  bool fIsExternalFunction = false;

  int  fNParams {}; // how many parameters for this gaincal
  TF1 fFitFunction;

  static int const MAXCHNS =   1;
  static int const MAXROWS =  80;
  static int const MAXCOLS =  52;
  static int const MAXROCS =   6;

  static int const NCHNS =   1;

  int const NROCS = 6;

  // Switched from float GC[NCHNS][NROCS][NCOLS][NROWS][6]
  std::vector<std::vector<std::vector<std::vector<std::vector<float > > > > > GC;
  std::vector<std::pair<float, float>> VC;  // VCal Calibration


  // Map for hardware locations by fed channel
  std::map<int, int> fHardwareMap;

};















#endif
