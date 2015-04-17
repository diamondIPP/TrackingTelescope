#include "PSIRootFileReader.h"

#include <iostream>
#include <string>
#include <stdint.h>
#include <stdlib.h>  



PSIRootFileReader::PSIRootFileReader (std::string const InFileName,
				      std::string const CalibrationList,
				      std::string const AlignmentFileName,
				      int const nrocs,
				      bool const useGainInterpolator
				      ) : PSIFileReader(CalibrationList,
							AlignmentFileName,
							nrocs,
							useGainInterpolator)
{
  
  fFileName = InFileName;
  if (!OpenFile()) {
    std::cerr << "ERROR: cannot open input file: " << fFileName << std::endl;
    throw;
  }
}


PSIRootFileReader::~PSIRootFileReader ()
{ 
  Clear(); 
  delete fTree;
  
  // Crashed when uncommented. Live with the memleak for now
  //delete fRootFile;
}


bool PSIRootFileReader::OpenFile ()
{
  
  fRootFile = new TFile(fFileName.c_str(), "READ");
  
  if (!fRootFile->IsOpen())
    return false;

  fTree = (TTree*)fRootFile->Get("tree");
  
  fAtEntry = 0;
  fNEntries = fTree->GetEntries();

  if (!(fNEntries>0))
    return false;

  // Set Branch Addresses
  fTree->SetBranchAddress("event_number",&f_event_number);
  fTree->SetBranchAddress("time",&f_time);

  fTree->SetBranchAddress("plane", &f_plane);
  fTree->SetBranchAddress("col", &f_col);
  fTree->SetBranchAddress("row", &f_row);
  fTree->SetBranchAddress("adc", &f_adc);
  fTree->SetBranchAddress("charge", &f_charge);

  return true;
}

void PSIRootFileReader::ResetFile ()
{  
  delete fTree;
  delete fRootFile;  
  OpenFile();
}

int PSIRootFileReader::GetNextEvent (){

  Clear(); 

  if (fAtEntry == fNEntries)
    return -1;

  fTree->GetEntry(fAtEntry);
  fAtEntry++;

  for (int i_hit = 0; i_hit != f_plane->size(); i_hit++){

    int roc = (*f_plane)[i_hit];
    int col = (*f_col)[i_hit];
    int row = (*f_row)[i_hit];
    int adc = (*f_adc)[i_hit];

    if (!IsPixelMasked( 1*100000 + roc*10000 + col*100 + row)){
      PLTHit* Hit = new PLTHit(1, roc, col, row, adc);
      
      // TODO: Either use calibration or read charge from root file
      //if (fUseGainInterpolator)
      //  fGainInterpolator.SetCharge(*Hit);
      //else
      //  fGainCal.SetCharge(*Hit);	
      
      // Dummy calibration for now
      Hit->SetCharge((*f_adc)[i_hit] * 65);
      
      fAlignment.AlignHit(*Hit);
      fHits.push_back(Hit);
      fPlaneMap[Hit->ROC()].AddHit(Hit);
    }
  }

  
  // Loop over all planes and clusterize each one, then add each plane to the correct telescope (by channel number
  for (std::map< int, PLTPlane>::iterator it = fPlaneMap.begin(); it != fPlaneMap.end(); ++it) {
    it->second.Clusterize(PLTPlane::kClustering_AllTouching, PLTPlane::kFiducialRegion_All);
    AddPlane( &(it->second) );
  }

  
  // If we are doing single plane-efficiencies:
  // Just send all events to the tracking and sort it out there
  if (DoingSinglePlaneEfficiency()){
    RunTracking( *((PLTTelescope*) this));
  }
  // Otherwise require exactly one hit per plane
  else {
    // 6 AND 3F FOR 6PLANES: TODO -> make dynamic!!!
    if (NClusters() == 4 && HitPlaneBits() == 0xf) {
      RunTracking( *((PLTTelescope*) this));
    }
  }
  
  return 0;

}
