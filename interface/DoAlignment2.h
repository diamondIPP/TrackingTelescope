#ifndef DoAlignment2_h
#define DoAlignment2_h

#include "PSIBinaryFileReader.h"
#include "PSIRootFileReader.h"
#include "GetNames.h"
#include "TestPlaneEfficiencySilicon.h"
#include "PLTPlane.h"
#include "PLTAlignment.h"

int DoAlignment2 (std::string, TFile*, TString, int);
void print_progress2(uint32_t, uint32_t);

#endif // DoAlignment2_h
