#ifndef TestBinaryFileReader_h
#define TestBinaryFileReader_h

#include <iostream>
#include <vector>
#include <string>
#include <ctime>

#include "TLegend.h"
#include "TLegendEntry.h"
#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH3F.h"
#include "TProfile2D.h"
#include "TParameter.h"
#include "TTree.h"
#include "TText.h"
#include "TPad.h"
#include "TROOT.h"

#include "PSIRootFileReader.h"
#include "PSIBinaryFileReader.h"
#include "GetNames.h"
#include "TestPlaneEfficiencySilicon.h"

int TestPSIBinaryFileReader (std::string const InFileName, TFile * out_f, TString const RunNumber,
                             int telescopeID);

void TestPlaneEfficiency (std::string const InFileName,
                          TFile * out_f,
                          TString const RunNumber,
                          int plane_under_test,
                          int n_events,
                          int telescopeID);

void WriteHTML (TString const OutDir, TString const CalFile, int telescopeID);


#endif //TestBinaryFileReader_h
