#ifndef GUARD_PSIRootFileReader_h
#define GUARD_PSIRootFileReader_h

#include "PSIFileReader.h"

#include "TFile.h"
#include "TTree.h"

class PSIRootFileReader : public PSIFileReader
{
  public:
    PSIRootFileReader (std::string const, 
		       std::string const, 
		       std::string const,
		       int const,
		       bool const,
		       bool const);
    ~PSIRootFileReader ();

    bool OpenFile ();
    void ResetFile ();
    int GetNextEvent ();

  private:
    std::string fFileName;
    TFile * fRootFile;
    TTree * fTree;

    //  Current entry and total number of entries
    int fAtEntry;
    int fNEntries;
  
    // Scalar Branches     
    int f_event_number;
    float f_time;

    // Vector Branches     
    std::vector<int> * f_plane;
    std::vector<int> * f_col;
    std::vector<int> * f_row;
    std::vector<int> * f_adc;
    std::vector<int> * f_charge;  
};

#endif
