#include "FileWriterTracking.h"
#include "GetNames.h"

using namespace std;

/** ============================
 CONSTRUCTOR
 =================================*/
FileWriterTracking::FileWriterTracking(string InFileName, uint8_t telescopeID, PSIFileReader * FR):
    nRoc(GetNumberOfROCS(telescopeID)) {

    NewFileName = getFileName(InFileName);
    br_charge_all.resize(nRoc);
    intree = ((PSIRootFileReader*) FR)->fTree;

}
FileWriterTracking::~FileWriterTracking(){ }

/** ============================
 AUXILIARY FUNCTIONS
 =================================*/
string FileWriterTracking::getFileName(string InFileName){
    stringstream ss(InFileName);
    while (getline(ss, NewFileName, '/')){}
    NewFileName.insert(int(NewFileName.length()-5), "_withTracks");
    return NewFileName;
}
