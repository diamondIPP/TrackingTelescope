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
    PSIRootFileReader(  std::string const,
                        std::string const,
                        std::string const,
                        int const,
                        bool const,
                        bool const,
                        bool const onlyAlign = false,
                        uint8_t const TelescopeID=0,
                        bool TrackOnlyTelescope=false);
    ~PSIRootFileReader ();

    bool OpenFile () override;
    void ResetFile () override;
    int GetNextEvent () override;
    void CloseFile() override;
    void ClearVectors();
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
    std::vector<uint16_t> * f_plane;
    std::vector<uint16_t> * f_col;
    std::vector<uint16_t> * f_row;
    std::vector<int16_t> * f_adc;
    std::vector<uint32_t> * f_charge;
    std::vector<float> * f_signal;
};

#endif
