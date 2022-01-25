#ifndef GUARD_PSIRootFileReader_h
#define GUARD_PSIRootFileReader_h

#include "PSIFileReader.h"

#include "TFile.h"
#include "TTree.h"
#include "TMacro.h"
#include "GetNames.h"
#include <iomanip>

class PSIRootFileReader : public PSIFileReader
{
  public:
    PSIRootFileReader(std::string in_file_name, bool only_align, bool track_only_telescope);
    ~PSIRootFileReader () override;

    bool OpenFile () override;
    void ResetFile () override;
    int GetNextEvent () override;
    void CloseFile() override;
    unsigned GetEntries() override { return fTree->GetEntries(); }

    // Make tree accessible
    TTree * fTree;
    TMacro * fMacro;
    TFile * fRootFile;

  private:
    std::string fFileName;

    const bool fOnlyAlign;

    //  Current entry and total number of entries
    int fAtEntry;
    int fNEntries;

    // Scalar Branches
    int32_t f_event_number;
    double f_time;

    // Vector Branches
    UShort_t f_n_hits;
    uint8_t f_plane[UINT8_MAX + 1] {};
    uint8_t f_col[UINT8_MAX + 1] {};
    uint8_t f_row[UINT8_MAX + 1] {};
    int16_t f_adc[UINT8_MAX + 1] {};
    float f_charge[UINT8_MAX + 1] {};
    float f_signal[UINT8_MAX + 1] {};
};

#endif
