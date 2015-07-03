#ifndef TestPlaneEfficiencySilicon_h
#define TestPlaneEfficiencySilicon_h

#include <string>
#include <vector>

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

#include "PSIRootFileReader.h"
#include "PSIBinaryFileReader.h"
#include "GetNames.h"

#define DEBUG false

int TestPlaneEfficiencySilicon (std::string const InFileName, TFile * out_f,
                                 TString const RunNumber, int telescopeID);


#endif // TestPlaneEfficiencySilicon_h
