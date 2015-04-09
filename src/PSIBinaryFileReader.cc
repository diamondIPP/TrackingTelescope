#include "PSIBinaryFileReader.h"

#include <iostream>
#include <string>
#include <stdint.h>
#include <stdlib.h>  

#include "TGraph.h"
#include "TString.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include "TMarker.h"
#include "TLine.h"

PSIBinaryFileReader::PSIBinaryFileReader (std::string const InFileName,
                                          std::string const CalibrationList,
                                          std::string const AlignmentFileName,
					  int const nrocs,
					  bool const useGainInterpolator
					  ) : NMAXROCS(nrocs), 
					      fGainCal(nrocs),
					      fUseGainInterpolator(useGainInterpolator){
  // constructor
  fEOF = 0;
  fBinaryFileName = InFileName;
  if (!OpenFile()) {
    std::cerr << "ERROR: cannot open input file: " << InFileName << std::endl;
    throw;
  }


  // Initialize fLevelsROC with zeros
  for (int i_roc=0; i_roc != NMAXROCS; i_roc++){    
    std::vector<float> tmp;
    for (int i_lev=0; i_lev != 6; i_lev++){
      tmp.push_back(0.);
    }
    fLevelsROC.push_back(tmp);    
  }

  // Initialize fCalibrationFile and fRawCalibrationFile with empty
  // strings
  for (int i_roc=0; i_roc != NMAXROCS; i_roc++){    
    fCalibrationFile.push_back("");
    fRawCalibrationFile.push_back("");
  }
  

  std::ifstream fCL(CalibrationList.c_str());
  if (!fCL.is_open()) {
      std::cerr << "ERROR: cannot open calibration list file: " << CalibrationList << std::endl;
      throw;
  }

  fCL >> fBaseCalDir;
  
  for (int i_roc=0; i_roc != NMAXROCS; i_roc++)
    fCL >> fCalibrationFile[i_roc];    

  for (int i_roc=0; i_roc != NMAXROCS; i_roc++)
    fCL >> fRawCalibrationFile[i_roc];

  for (int i_roc=0; i_roc != NMAXROCS; i_roc++)
    fGainCal.ReadGainCalFile(fBaseCalDir + "/" + fCalibrationFile[i_roc]);

  for (int i_roc=0; i_roc != NMAXROCS; i_roc++)
    fGainInterpolator.ReadFile(fBaseCalDir + "/" + fRawCalibrationFile[i_roc], i_roc);

  fAlignment.ReadAlignmentFile(AlignmentFileName);
  SetTrackingAlignment(&fAlignment);
  
  // TODO: Make sure this works correctly also for <6 planes!!
  SetTrackingAlgorithm(PLTTracking::kTrackingAlgorithm_6PlanesHit);

}


PSIBinaryFileReader::~PSIBinaryFileReader ()
{
  // Destructor
  Clear();
 
}


bool PSIBinaryFileReader::OpenFile ()
{

  fEOF = false;
  fInputBinaryFile.clear();
  fInputBinaryFile.open(fBinaryFileName.c_str(), std::ios::in | std::ios::binary);
  if (fInputBinaryFile.is_open()) {
    return true;
  }

  return false;
}

void PSIBinaryFileReader::ResetFile ()
{
  // Reset the file
  std::cout << "Reset the file!" << std::endl;

  fInputBinaryFile.clear() ;
  fInputBinaryFile.seekg(0, fInputBinaryFile.beg) ;
  fEOF = false;
 
//  sleep(10);
//
//  OpenFile();  
}



unsigned short PSIBinaryFileReader::readBinaryWordFromFile()
{

  if (fInputBinaryFile.eof()) fEOF = 1;

  unsigned char a = fInputBinaryFile.get();
  unsigned char b = fInputBinaryFile.get();
  unsigned short word  =  (b << 8) | a;

//    cout << Form("readBinaryWord: %02x %02x word: %04x ", a, b, word) << endl;

  return word;
}

// ----------------------------------------------------------------------
int PSIBinaryFileReader::nextBinaryHeader()
{

  int header(-1);

  for (int i = 0; i < MAXNDATA; ++i)
  {
    fBuffer[i] = 0;
  }

  fBufferSize = 0;
  unsigned short word(0);

  while (1)
  {
    word = readBinaryWordFromFile();

    if (fEOF) break;

    if (word == 0x8000)
    {
      //      cout << "==> end of file " << endl;
      header = 0;
      break;
    }

    if (word == 0x8001)
    {
      header = 1;
      //      cout << " --> data " <<endl;
    }

    if (word == 0x8004)
    {
      header = 4;
      //      cout << " --> trig " << endl;
    }

    if (word == 0x8008)
    {
      header = 8;
      //      cout << " --> reset " << endl;
    }

    if (word == 0x8080)
    {
      header = 80;
      //      cout << " --> overflow " << endl;
    }

    if (header > -1)
    {
//       cout << " --> Found next header " << header << endl;
      break;
    }

    if (0) std::cout << Form("Adding %04x", word) << " at " << fBufferSize << std::endl;
    fBuffer[fBufferSize] = word;
    ++fBufferSize;
    if (fBufferSize >= MAXNDATA) {
      std::cerr << "ERROR: fBufferSize >= MAXNDATA: " << fBufferSize << std::endl;
      exit(1);
    }
  }

  fHeader     = fNextHeader;
  fNextHeader = header;

  return header;
}


// ----------------------------------------------------------------------
int PSIBinaryFileReader::decodeBinaryData()
{
  int j(0);

  if (fHeader > 0)  {
    unsigned short t0 = fBuffer[0];
    unsigned short t1 = fBuffer[1];
    unsigned short t2 = fBuffer[2];
    fUpperTime = t0;
    fLowerTime = (t1 << 16) | t2;
    fTime = (((long long int) fUpperTime)<<32) + ((long long int)fLowerTime);
    //    cout << Form(" Event at time  %04x/%08x with Header %d", fUpperTime, fLowerTime, fHeader) << endl;
    //    cout << Form(" Event at time  %04d/%08d with Header %d", fUpperTime, fLowerTime, fHeader) << endl;
  } else  {
    //    cout << "No valid header, skipping this event" << endl;
    return -1;
  }

  int value(0);
  for (int i = 3; i < fBufferSize; ++i)
  {
    value = fBuffer[i] & 0x0fff;
    if (value & 0x0800) value -= 4096;
    fData[i-3] = value;
    ++j;
  }

  fBufferSize -=3;

  if (0)
  {
    for (int i = 0; i < fBufferSize; ++i)
    {
      std::cout << Form(" %04x ", fData[i]);
    }
    std::cout << std::endl;

    for (int i = 0; i < fBufferSize; ++i)
    {
      std::cout << Form(" %6i ", fData[i]);
    }
    std::cout << std::endl;
  }

  return j;

}




size_t PSIBinaryFileReader::NHits ()
{
  return fHits.size();
}


PLTHit* PSIBinaryFileReader::Hit (size_t const i)
{
  if (i < fHits.size()) {
    return fHits[i];
  }
  std::cerr << "ERROR: PSIBinaryFileReader::Hit asking for a hit outside of range." << std::endl;
  throw;
}


void PSIBinaryFileReader::Clear()
{
  for (size_t i = 0; i != fHits.size(); ++i) {
    if (fHits[i]) {
      delete fHits[i];
    }
  }
  fHits.clear();

  fPlaneMap.clear();
  fPlanes.clear();
  fTracks.clear();

  return;
}


int PSIBinaryFileReader::GetNextEvent ()
{
  Clear();

  // Simple cheat to always have 6 planes
  // (well, 4 for now)
  for (int i = 0; i != 4; ++i) {
    fPlaneMap[i].SetROC(i);;
  }

  while (nextBinaryHeader() >= 0) {
    decodeBinaryData();
    if (fBufferSize <= 0) {
      continue;
    } else {
      // decode waveform
      DecodeHits();
      return fBufferSize;
    }
  }

  // no more events, return -1
  return -1;
}




int PSIBinaryFileReader::CalculateLevels (int const NMaxEvents,TString const OutDir)
{
  // Vector for ROC level histograms
  std::vector<TH1F*> hROCLevels;
  bool CreateLevelHists = true;


  // Hist for TBM Levels
  TH1F hTBMLevels("LevelsTBM", "LevelsTBM", 100, -1000, 1000);

  //for (int ievent = 0; NMaxEvents > 0 ? ievent != NMaxEvents : !fEOF; ++ievent) {
  for (int ievent = 0; !fEOF; ++ievent) {
    while (nextBinaryHeader() >= 0) {
      decodeBinaryData();
      if (fBufferSize <= 0) {
        continue;
      }
        break;
    }
      if (fBufferSize <= 0) {
        continue;
      }

    // Now I should have a waveform in fData of length fBufferSize
    int UBCount = 0;
    std::vector<int> UBPosition;

    for (int i = 0; i != fBufferSize; ++i) {
      if (fData[i] < UBLevel) {
        ++UBCount;
        UBPosition.push_back(i);
      }
    }

    int const NROCs = UBCount - 5;
    //std::cout << "fBufferSize: " << fBufferSize << std::endl;
    //std::cout << "NROCs: " << NROCs << std::endl;
    if (NROCs != NMAXROCS) {
      std::cerr << "WARNING: NROCs != NMAXROCS in levels calculation.  Skipping event" << std::endl;
      continue;
    }
    if (NROCs <= 0) {
      std::cerr << "ERROR: bad event with NROCs <= 0" << std::endl;
      continue;
    }

    if (CreateLevelHists) {
      CreateLevelHists = false;
      hROCLevels.resize(NROCs);
      for (int ihist = 0; ihist != NROCs; ++ihist) {
        TString const Name = TString::Format("Levels_ROC%i", ihist);
        std::cout << "Creating level hist: " << Name << std::endl;
        hROCLevels[ihist] = new TH1F(Name.Data(), Name.Data(), 200, -1000, 1000);
      }
    }

    if (NROCs > (int) hROCLevels.size()) {
      std::cerr << "ERROR: For some reason this event has more rocs than before.  Skipping it" << std::endl;
      continue;
    }

    hTBMLevels.Fill( fData[UBPosition[0] + 4] );
    hTBMLevels.Fill( fData[UBPosition[0] + 5] );
    hTBMLevels.Fill( fData[UBPosition[0] + 6] );
    hTBMLevels.Fill( fData[UBPosition[0] + 7] );


    std::vector<int> UBPositionROC(NROCs, -1);
    std::vector<int> NHitsROC(NROCs, 0);
    for (int iroc = 0; iroc != NROCs; ++iroc) {
      UBPositionROC[iroc] = UBPosition[3+iroc];
      NHitsROC[iroc] = (UBPosition[3+iroc+1] - UBPosition[3+iroc] - 3) / 6;
      //if (NHitsROC[iroc] > 0) {
      //  printf("ievent: %5i  roc %2i  NHits: %5i\n", ievent, iroc, NHitsROC[iroc]);
      //}

      for (int ihit = 0; ihit != NHitsROC[iroc]; ++ihit) {
        hROCLevels[iroc]->Fill( fData[ UBPositionROC[iroc] +  + 2 + 1 + ihit * 6]);
        hROCLevels[iroc]->Fill( fData[ UBPositionROC[iroc] +  + 2 + 2 + ihit * 6]);
        hROCLevels[iroc]->Fill( fData[ UBPositionROC[iroc] +  + 2 + 3 + ihit * 6]);
        hROCLevels[iroc]->Fill( fData[ UBPositionROC[iroc] +  + 2 + 4 + ihit * 6]);
        hROCLevels[iroc]->Fill( fData[ UBPositionROC[iroc] +  + 2 + 5 + ihit * 6]);
      }

    }

  }


  for (size_t iroc = 0; iroc != hROCLevels.size(); ++iroc) {
    TSpectrum Spectrum(20);
    Spectrum.SetAverageWindow(30);//probably does nothing
    int const NPeaks = Spectrum.Search(hROCLevels[iroc]);
    float* Peaks = Spectrum.GetPositionX();
    std::sort(Peaks, Peaks + NPeaks);

    // Workaround for:
    //  TMarker pPoint[NPeaks];
    std::vector< TMarker> pPoint(NPeaks);
    for (int i = 0; i < NPeaks; ++i) {
      pPoint[i].SetX(Peaks[i]);
      pPoint[i].SetY(hROCLevels[iroc]->GetBinContent(hROCLevels[iroc]->FindBin(Peaks[i])));
      pPoint[i].SetMarkerStyle(3);
    }

    float const hHistMaxY = hROCLevels[iroc]->GetMaximum();

    // Workaround for
    // TLine lLine[NPeaks];
    std::vector< TLine > lLine(NPeaks);

    printf("Levels calculated for ROC %i:\n", (int) iroc);
    for (int i = 0; i < NPeaks - 1; ++i) {
      float xp = Peaks[i];
      float yp = Peaks[i + 1];
      xp = xp + (yp - xp) / 2.0;

      printf("  Threshold %d value %f\n", i, xp);

      if (i < 6) {
        fLevelsROC[iroc][i] = xp;
      }

      lLine[i].SetLineColor(2);
      lLine[i].SetX1(xp);  lLine[i].SetX2(xp);
      lLine[i].SetY1(1);   lLine[i].SetY2(hHistMaxY);
    }




    TCanvas Can;
    Can.cd();
    hROCLevels[iroc]->Draw("hist");
    for (int i = 0; i < NPeaks; ++i) {
      pPoint[i].Draw("same");
      lLine[i].Draw("same");
    }
    Can.SaveAs(OutDir + "/" + TString(hROCLevels[iroc]->GetName()) + ".gif");
  }

  TCanvas Can;
  Can.cd();
  hTBMLevels.Draw("hist");
  Can.SaveAs(OutDir + "LevelsTBM.gif");

  ResetFile();

  return 0;
}



void PSIBinaryFileReader::DecodeHits ()
{
  int UBCount = 0;
  std::vector<int> UBPosition;

  for (int i = 0; i != fBufferSize; ++i) {
    if (fData[i] < UBLevel) {
      ++UBCount;
      UBPosition.push_back(i);
    }
  }

  int const NROCs = UBCount - 5;
  //std::cout << "fBufferSize: " << fBufferSize << std::endl;
  //std::cout << "NROCs: " << NROCs << std::endl;
  //static int iBadData = 0;
  //static int iGoodData = 0;
  if (NROCs > NMAXROCS) {
    //DrawWaveform(TString::Format("BadWave_%i.gif", iBadData++));
    std::cerr << "WARNING: NROCs > NMAXROCS.  Too many UBs.  Skipping this event" << std::endl;
    return;
  } else if (NROCs != NMAXROCS) {
    //DrawWaveform(TString::Format("BadWave_%i.gif", iBadData++));
    std::cerr << "WARNING: NROCs != NMAXROCS.  Skipping this event" << std::endl;
    return;
  }
  //DrawWaveform(TString::Format("GoodWave_%i.gif", iGoodData++));
  if (NROCs <= 0) {
    std::cerr << "WARNING: bad event with NROCs <= 0" << std::endl;
    return;
  }






  std::vector<int> UBPositionROC(NROCs, -1);
  std::vector<int> NHitsROC(NROCs, 0);
  for (int iroc = 0; iroc != NROCs; ++iroc) {
    UBPositionROC[iroc] = UBPosition[3+iroc];
    NHitsROC[iroc] = (UBPosition[3+iroc+1] - UBPosition[3+iroc] - 3) / 6;
    //printf("roc %2i  NHits: %5i\n", iroc, NHitsROC[iroc]);

    for (int ihit = 0; ihit != NHitsROC[iroc]; ++ihit) {
      std::pair<int, int> colrow = fill_pixel_info(fData, UBPosition[3 + iroc] + 2 + 0 + ihit * 6, iroc);

      // Ignore masked pixels
      // Important: Assume Channel==1 !!!!
      if (!IsPixelMasked( 1*100000 + iroc*10000 + colrow.first*100 + colrow.second)){
        //printf("Hit iroc %2i  col %2i  row %2i  PH: %4i\n", iroc, colrow.first, colrow.second, fData[ UBPosition[3 + iroc] + 2 + 6 + ihit * 6 ]);
        PLTHit* Hit = new PLTHit(1, iroc, colrow.first, colrow.second, fData[ UBPosition[3 + iroc] + 2 + 6 + ihit * 6 ]);

	if (fUseGainInterpolator)
	  fGainInterpolator.SetCharge(*Hit);
	else
	  fGainCal.SetCharge(*Hit);	

        fAlignment.AlignHit(*Hit);
        fHits.push_back(Hit);
        fPlaneMap[Hit->ROC()].AddHit(Hit);
      }
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

  return;
}





int PSIBinaryFileReader::LevelInfo (int const Value, int const iroc)
{
  //if (Value <= 0) {
    //std::cout << "Something is wrong" << std::endl;
    //return -1;
  //}
       if (                               Value <= fLevelsROC[iroc][0]) return 0;
  else if (Value > fLevelsROC[iroc][0] && Value <= fLevelsROC[iroc][1]) return 1;
  else if (Value > fLevelsROC[iroc][1] && Value <= fLevelsROC[iroc][2]) return 2;
  else if (Value > fLevelsROC[iroc][2] && Value <= fLevelsROC[iroc][3]) return 3;
  else if (Value > fLevelsROC[iroc][3] && Value <= fLevelsROC[iroc][4]) return 4;
  else if (Value > fLevelsROC[iroc][4]) return 5;

  return -1;
}



std::pair<int, int> PSIBinaryFileReader::fill_pixel_info(int* evt , int ctr, int iroc)
{
  int finalcol = -1;
  int finalrow = -1;
  int c1 = evt[ctr + 1];
  int c0 = evt[ctr + 2];
  int r2 = evt[ctr + 3];
  int r1 = evt[ctr + 4];
  int r0 = evt[ctr + 5];

  int trancol=(LevelInfo(c1, iroc))*6 + (LevelInfo(c0, iroc));
  int tranrow=(LevelInfo(r2, iroc))*36 + (LevelInfo(r1, iroc))*6 + (LevelInfo(r0, iroc));
  if(tranrow%2 == 0)
  {
    finalrow=79 - (tranrow - 2)/2;
    finalcol=trancol * 2;
  }
  else
  {
    finalrow=79 - (tranrow - 3)/2;
    finalcol=trancol * 2 + 1;
  }



  return std::make_pair<int, int>(finalcol, finalrow);
}


void PSIBinaryFileReader::AddToPixelMask(int ch, int roc, int col, int row)
{
  fPixelMask.insert( ch*100000 + roc*10000 + col*100 + row );
}


void PSIBinaryFileReader::ReadPixelMask (std::string const InFileName)
{
  std::cout << "PLTBinaryFileReader::ReadPixelMask reading file: " << InFileName << std::endl;

  std::ifstream InFile(InFileName.c_str());
  if (!InFile.is_open()) {
    std::cerr << "ERROR: cannot open PixelMask file: " << InFileName << std::endl;
    throw;
  }

  // Loop over header lines in the input data file
  for (std::string line; std::getline(InFile, line); ) {
    if (line == "") {
      break;
    }
  }

  for (std::string line; std::getline(InFile, line); ) {
    int ch=0, roc=0, col=0, row=0;
    std::istringstream linestream;

    linestream.str(line);
    linestream >> ch >> roc >> col >> row;

    //std::cout << "Masking: " << ch << " " << roc << " " << col << " " << row << std::endl;
    fPixelMask.insert( ch*100000 + roc*10000 + col*100 + row );
  }

  return;
}


bool PSIBinaryFileReader::IsPixelMasked (int const ChannelPixel)
{
  if (fPixelMask.count(ChannelPixel)) {
    return true;
  }
  return false;
}



bool PSIBinaryFileReader::ReadAddressesFromFile (std::string const InFileName)
{
  // Open file
  std::ifstream f(InFileName.c_str());
  if (!f.is_open()) {
    std::cerr << "ERROR: PSIBinaryFileReader::ReadAddressesFromFile cannot open file: " << InFileName << std::endl;
    throw;
  }

  std::cout << "PSIBinaryFileReader::ReadAddressesFromFile: " << InFileName << std::endl;

  std::string line;
  std::getline(f, line);
  std::getline(f, line);
  std::getline(f, line);
  std::getline(f, line);
  std::getline(f, line);
  std::getline(f, line);
  std::getline(f, line);
  std::getline(f, line);
  std::getline(f, line);

  std::string ROC;
  int unused;
  while(f >> ROC >> unused >> unused >> unused) {
    int const iroc = atoi( ROC.substr(3, ROC.length()-1).c_str());
    f >> fLevelsROC[iroc][0]
      >> fLevelsROC[iroc][1]
      >> fLevelsROC[iroc][2]
      >> fLevelsROC[iroc][3]
      >> fLevelsROC[iroc][4]
      >> fLevelsROC[iroc][5]
      >> unused;

    printf("Set Levels ROC %2i :  ", iroc);
    for (int i = 0; i != 6; ++i) {
      printf("  %4i", (int) fLevelsROC[iroc][i]);
    }
    printf("\n");
  }

  return true;
}





void PSIBinaryFileReader::DrawTracksAndHits (std::string const Name)
{
  int const NH = NHits();
  int const NC = NClusters();
  int const NT = NTracks();

  float X[NH];
  float Y[NH];
  float Z[NH];

  float CX[NC];
  float CY[NC];
  float CZ[NC];

  // Workaround for:
  //  TLine Line[3][NT];
  // which does not work in CLANG
  std::vector< std::vector< TLine > > Line;
  for (int i=0; i<3; i++){
    std::vector <TLine> tmp(NT);
    Line.push_back( tmp );
  }

  TH2F* HistCharge[6];
  for (int i = 0; i != 6; ++i) {
    TString Name;
    Name.Form("ChargeMap_ROC%i", i);
    HistCharge[i] = new TH2F(Name, Name, PLTU::NCOL, PLTU::FIRSTCOL, PLTU::LASTCOL+1, PLTU::NROW, PLTU::FIRSTROW, PLTU::LASTROW);
    HistCharge[i]->GetZaxis()->SetRangeUser(0, 50000);
  }

  TH2F* HistChargeUnclustered[6];
  for (int i = 0; i != 6; ++i) {
    TString Name;
    Name.Form("ChargeMapUnclustered_ROC%i", i);
    HistChargeUnclustered[i] = new TH2F(Name, Name, PLTU::NCOL, PLTU::FIRSTCOL, PLTU::LASTCOL+1, PLTU::NROW, PLTU::FIRSTROW, PLTU::LASTROW);
    HistChargeUnclustered[i]->GetZaxis()->SetRangeUser(0, 50000);
  }



  int j = 0;
  for (size_t ip = 0; ip != NPlanes(); ++ip) {
    PLTPlane* P = Plane(ip);
    for (size_t ih = 0; ih != P->NHits(); ++ih) {
      PLTHit* H = P->Hit(ih);
      X[j] = H->TX();
      Y[j] = H->TY();
      Z[j] = H->TZ();
      ++j;

      HistCharge[ip]->SetBinContent(H->Column() + 1 - PLTU::FIRSTCOL, H->Row() + 1 - PLTU::FIRSTROW, H->Charge());
    }
    for (size_t ih = 0; ih != P->NUnclusteredHits(); ++ih) {
      PLTHit* H = P->UnclusteredHit(ih);
      HistChargeUnclustered[ip]->SetBinContent(H->Column() + 1 - PLTU::FIRSTCOL, H->Row() + 1 - PLTU::FIRSTROW, H->Charge());
    }



  }
  int jc = 0;
  for (size_t ip = 0; ip != NPlanes(); ++ip) {
    PLTPlane* P = Plane(ip);
    for (size_t ic = 0; ic != P->NClusters(); ++ic) {
      PLTCluster* C = P->Cluster(ic);
      CX[jc] = C->TX();
      CY[jc] = C->TY();
      CZ[jc] = C->TZ();
      ++jc;
    }
  }

  std::vector<PLTHit*> UsedHits;

  for (int i = 0; i != NT; ++i) {
    PLTTrack* T = fTracks[i];

    // XZ
    Line[0][i].SetX1(0);
    Line[0][i].SetX2(10.16);
    Line[0][i].SetY1(T->TX(0));
    Line[0][i].SetY2(T->TX(10.16));
    Line[0][i].SetLineColor(i+1);

    // YZ
    Line[1][i].SetX1(0);
    Line[1][i].SetX2(10.16);
    Line[1][i].SetY1(T->TY(0));
    Line[1][i].SetY2(T->TY(10.16));
    Line[1][i].SetLineColor(i+1);

    // XY
    Line[2][i].SetX1(T->TX(0));
    Line[2][i].SetX2(T->TX(10.16));
    Line[2][i].SetY1(T->TY(0));
    Line[2][i].SetY2(T->TY(10.16));
    Line[2][i].SetLineColor(i+1);

    //printf("XY 0 7: %9.3f %9.3f   %9.3f %9.3f\n", T->TX(0), T->TY(0), T->TX(7), T->TY(7));
  }

  //for (int i = 0; i != 3; ++i) {
  //  for (int j = 0; j != NT; ++j) {
  //    Line[i][j].SetLineColor(4);
  //  }
  //}

  TCanvas C("TelescopeTrack", "TelescopeTrack", 800, 800);;
  C.Divide(3, 3);

  C.cd(1);
  TGraph gXZ(NC, CZ, CX);
  gXZ.SetTitle("");
  gXZ.GetXaxis()->SetTitle("Z (cm)");
  gXZ.GetYaxis()->SetTitle("X (cm)");
  gXZ.GetXaxis()->SetTitleSize(0.06);
  gXZ.GetYaxis()->SetTitleSize(0.08);
  gXZ.GetXaxis()->SetTitleOffset(0.7);
  gXZ.GetYaxis()->SetTitleOffset(0.5);
  gXZ.SetMarkerColor(40);
  gXZ.GetXaxis()->SetLimits(-0.5, 11.0);
  gXZ.SetMinimum(-0.5);
  gXZ.SetMaximum( 0.5);
  if (NC) {
    gXZ.Draw("A*");
  }
  for (int i = 0; i != NT; ++i) {
    Line[0][i].Draw();
  }

  C.cd(4);
  TGraph gYZ(NC, CZ, CY);
  gYZ.SetTitle("");
  gYZ.GetXaxis()->SetTitle("Z (cm)");
  gYZ.GetYaxis()->SetTitle("Y (cm)");
  gYZ.GetXaxis()->SetTitleSize(0.06);
  gYZ.GetYaxis()->SetTitleSize(0.08);
  gYZ.GetXaxis()->SetTitleOffset(0.7);
  gYZ.GetYaxis()->SetTitleOffset(0.5);
  gYZ.SetMarkerColor(40);
  gYZ.GetXaxis()->SetLimits(-0.5, 11.0);
  gYZ.SetMinimum(-0.5);
  gYZ.SetMaximum( 0.5);
  if (NC) {
    gYZ.Draw("A*");
  }
  for (int i = 0; i != NT; ++i) {
    Line[1][i].Draw();
  }

  //TVirtualPad* Pad = C.cd(3);
  //Pad->DrawFrame(-30, -30, 30, 30);
  C.cd(7);
  TGraph gXY(NC, CX, CY);
  gXY.SetTitle("");
  gXY.GetXaxis()->SetTitle("X (cm)");
  gXY.GetYaxis()->SetTitle("Y (cm)");
  gXY.GetXaxis()->SetTitleSize(0.06);
  gXY.GetYaxis()->SetTitleSize(0.08);
  gXY.GetXaxis()->SetTitleOffset(0.7);
  gXY.GetYaxis()->SetTitleOffset(0.5);
  gXY.SetMarkerColor(40);
  gXY.GetXaxis()->SetLimits(-0.5, 0.5);
  gXY.SetMinimum(-0.5);
  gXY.SetMaximum( 0.5);
  if (NC) {
    gXY.Draw("A*");
  }
  for (int i = 0; i != NT; ++i) {
    Line[2][i].Draw();
  }

  C.cd(2);
  HistCharge[0]->Draw("colz");
  C.cd(3);
  HistCharge[1]->Draw("colz");
  C.cd(5);
  HistCharge[2]->Draw("colz");
  C.cd(6);
  HistCharge[3]->Draw("colz");
  C.cd(8);
  HistCharge[4]->Draw("colz");
  C.cd(9);
  HistCharge[5]->Draw("colz");

  //C.cd(3);
  //HistChargeUnclustered[0]->Draw("colz");
  //C.cd(6);
  //HistChargeUnclustered[1]->Draw("colz");
  //C.cd(9);
  //HistChargeUnclustered[2]->Draw("colz");

  C.SaveAs(Name.c_str());

  for (int i = 0; i != 6; ++i) {
    delete HistCharge[i];
    delete HistChargeUnclustered[i];
  }

  return;
}




void PSIBinaryFileReader::DrawWaveform (TString const OutFileName)
{
  int X[fBufferSize];
  for (int i = 0; i != fBufferSize; ++i) {
    X[i] = i;
  }
  TGraph g(fBufferSize, X, fData);
  TCanvas c;
  c.cd();
  g.Draw("AC*");
  c.SaveAs(OutFileName);

  return;
}
