#ifndef PLTANALYSIS_H
#define PLTANALYSIS_H

#include "PSIRootFileReader.h"
#include "PSIBinaryFileReader.h"
#include "GetNames.h"
#include "TestPlaneEfficiencySilicon.h"
#include "RootItems.h"
#include "FileWriterTracking.h"

#define verbose 0

/** ============================
 DEFAULT ANALYSIS CLASS
 =================================*/
class PLTAnalysis
{

private:
    uint8_t const telescopeID;
    std::string const InFileName;
    TString const RunNumber;
    TFile * out_f;
    /** measure elapsed time */
    float now, now1, now2, loop, startProg, endProg, allProg, averTime, speed;
    /** times for counting */
    uint32_t const TimeWidth, StartTime;
    uint32_t ThisTime;
    uint16_t NGraphPoints;
    /** miscellaneous */
    uint32_t const PHThreshold;
    bool do_slope;
    PSIFileReader * FR;
    uint32_t nEntries;
    RootItems * Histos;
    FileWriterTracking * FW;

public:

    /** ============================
     CONSTRUCTOR // DECONSTRUCTOR
     =================================*/
    PLTAnalysis(std::string const inFileName, TFile * Out_f,  TString const runNumber, uint8_t const TelescopeID);
    ~PLTAnalysis();


    /** ============================
     EVENT LOOP
     =================================*/
    void EventLoop();


    /** ============================
     AFTER LOOP -> FINISH
     =================================*/
    void FinishAnalysis();


    /** ============================
     AUXILIARY FUNCTIONS
     =================================*/
    float getTime(float now, float & time);
    void SinglePlaneStudies();
    void InitFileReader();
    void PrintProcess(uint32_t);
    void MeasureSpeed(uint32_t);
    void WriteTrackingTree();
    void MakeAvgPH();
    void DrawTracks();
    void FillPHHistos(uint8_t, PLTCluster*);
    void FillOccupancyHiLo(PLTCluster*);
    void FillOccupancy(PLTPlane*);
    void FillOfflinePH(PLTTrack*, PLTCluster*);
};

void TestPlaneEfficiency (std::string const InFileName,
                          TFile * out_f,
                          TString const RunNumber,
                          int plane_under_test,
                          int n_events,
                          int telescopeID);

void WriteHTML (TString const OutDir, TString const CalFile, int telescopeID);
#endif // PLTANALYSIS_H