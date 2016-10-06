#include "PLTGainCal.h"

using namespace std;


PLTGainCal::PLTGainCal (): NROCS(6)
{
  ResetGC();
  fIsGood = false;
  fIsExternalFunction = false;
}

PLTGainCal::PLTGainCal (int nrocs,
			bool isExternalFunction ): NROCS(nrocs)
{
  ResetGC();
  fIsGood = false;
  fIsExternalFunction = isExternalFunction;
}

PLTGainCal::PLTGainCal (std::string const GainCalFileName, int const NParams): NROCS(6)
{
    ResetGC();
    fIsGood = false;
    fIsExternalFunction = false;
    fNParams = NParams;
    if (NParams == 5)
        ReadGainCalFile5(GainCalFileName);
    else if (NParams == 3)
        ReadGainCalFile3(GainCalFileName);
    else {
        std::cerr << "ERROR: I have no idea how many params you have" << std::endl;
        throw;
    }
}


PLTGainCal::~PLTGainCal ()
{
}


int PLTGainCal::RowIndex (int const i)
{
  return i - PLTU::FIRSTROW;
}


int PLTGainCal::ColIndex (int const i)
{
  return i - PLTU::FIRSTCOL;
}


int PLTGainCal::ChIndex (int const i)
{
  return i - 1;
}


int PLTGainCal::RocIndex (int const i)
{
  return i;
}


float PLTGainCal::GetCoef(int const i, int const ch, int const roc, int const col, int const row)
{
  // Get a coef, note roc number is 0, 1, 2
  int irow = RowIndex(row);
  int icol = ColIndex(col);
  int ich  = ChIndex(ch);
  int iroc = RocIndex(roc);
  if (irow < 0 || icol < 0 || ich < 0 || iroc < 0) {
    return -9999;
  }

  return GC[ich][iroc][icol][irow][i];
}


void PLTGainCal::SetCharge (PLTHit& Hit, uint8_t telescopeID)
{
  Hit.SetCharge( GetCharge(Hit.Channel(), telescopeID, Hit.ROC(), Hit.Column(), Hit.Row(), Hit.ADC()) );
  return;
}



float PLTGainCal::GetCharge(int const ch, int telescopeID, int const roc, int const col, int const row, int INadc)
{
    /** Get charge, note roc number is 0, 1, 2... */
    int const adc = INadc;
    if (ChIndex(ch)   >= MAXCHNS) { printf("ERROR: over MAXCHNS: %i\n", ch); };
    if (RowIndex(row) >= MAXROWS) { printf("ERROR: over MAXROWS: %i\n", row); };
    if (ColIndex(col) >= MAXCOLS) { printf("ERROR: over MAXCOLS: %i\n", col); };

    int16_t irow = RowIndex(row);
    int16_t icol = ColIndex(col);
    int16_t ich  = ChIndex(ch);
    int16_t iroc = RocIndex(roc);

    if (irow < 0 || icol < 0 || ich < 0 || iroc < 0) {
        return -9999;
    }

    float charge = -9999;

    /** external calibration */
    if (fIsExternalFunction) {
        for (int ipar = 0; ipar < fNParams; ++ipar) {
            fFitFunction.SetParameter(ipar, GC[ich][iroc][icol][irow][ipar]);
        }

        /** change calibration factor for digital planes*/
        if (UseDigitalCalibration(telescopeID)){
            if (iroc == 4 || iroc == 5){
                charge = 47 * fFitFunction.GetX(adc);
//                cout << "iroc: adc/charge: " << int(iroc) <<": " << adc << "/" << charge << endl;
            }
            else if (iroc == 6)         charge = 43.13 * fFitFunction.GetX(adc) + 333.0;
            else                        charge = 65. * fFitFunction.GetX(adc);

        }
        else                            charge = 65. * fFitFunction.GetX(adc);
    }
    /** old calibration */
    else {
        if (fNParams == 3) {
            charge = 65. * (float(adc * adc) * GC[ich][iroc][icol][irow][2] + float(adc) * GC[ich][iroc][icol][irow][1] + GC[ich][iroc][icol][irow][0]);

        }
        else if (fNParams == 5) {
            charge = 65. * (TMath::Power( (float) adc, 2) * GC[ich][iroc][icol][irow][0] + (float) adc * GC[ich][iroc][icol][irow][1] + GC[ich][iroc][icol][irow][2]
              + (GC[ich][iroc][icol][irow][4] != 0 ? TMath::Exp( (adc - GC[ich][iroc][icol][irow][3]) / GC[ich][iroc][icol][irow][4] ) : 0)
              );
        }
        else {
            std::cerr << "ERROR: PLTGainCal::GetCharge() I do not know of that number of fNParams: " << fNParams << std::endl;
            exit(1);
        }
    }

    /** Printing */
    if (PLTGainCal::DEBUGLEVEL) {
            printf("%2i %1i %2i %2i %4i %10.1f\n", ch, roc, col, row, adc, charge);
    }

    return charge;
}

void PLTGainCal::ReadGainCalFile (std::string const GainCalFileName, int roc)
{

  std::cout << "ReadGainCalFile" << std::endl;

  if (GainCalFileName == "") {
    fIsGood = false;
    return;
  }

  std::ifstream InFile(GainCalFileName.c_str());
//  std::cout << std::getline(InFile, line) << std::endl;
  if (!InFile.is_open()) {
    std::cerr << "ERROR: cannot open gaincal file: " << GainCalFileName << std::endl;
    throw;
  }
  TString CheckFirstLine;
  CheckFirstLine.ReadLine(InFile);

  if (CheckFirstLine.BeginsWith("Parameters of the vcal vs. pulse height fits")) {// DA: TODO else condition?
    std::cout << "PLTGainCal setting fIsExternalFunction" << std::endl;
    fIsExternalFunction = true;
  }

  // Loop over header lines in the input data file
  std::string line;
  for (std::string line; std::getline(InFile, line); ) {
    if (line == "") {
      break;
    }
  }




  //std::string line;
  std::getline(InFile, line);
  std::istringstream linestream;
  linestream.str(line);
  int i = 0;
  if (fIsExternalFunction) {
  } else {
   i -= 4;
  }
  for (float junk; linestream >> junk; ++i) {
  }
  InFile.close();
  fNParams = i;

  printf("PLTGainCal sees a parameter file with %i params\n", fNParams);

  if (fIsExternalFunction) {
    ReadGainCalFileExt(GainCalFileName, roc);
  } else {
    if (fNParams == 5) {
      ReadGainCalFile5(GainCalFileName);
    } else if (fNParams == 3) {
      ReadGainCalFile3(GainCalFileName);
    } else {
      std::cerr << "ERROR: I have no idea how many params you have" << std::endl;
      throw;
    }
  }

  return;
}

int PLTGainCal::GetHardwareID (int const Channel)
{
  return fHardwareMap[Channel];
}

void PLTGainCal::ReadGainCalFile5 (std::string const GainCalFileName)
{
  int ch, row, col, roc;
  int irow;
  int icol;
  int ich;

  std::cout << "Reading GainCal file with 5 parameters" << std::endl;

  ifstream f(GainCalFileName.c_str());
  if (!f) {
    std::cerr << "ERROR: cannot open file: " << GainCalFileName << std::endl;
    throw;
  }

  // Loop over header lines in the input data file
  for (std::string line; std::getline(f, line); ) {
    int mf, mfc, hub;
    if (line == "") {
      break;
    }
    std::istringstream ss;
    ss.str(line);
    ss >> mf >> mfc >> hub >> ch;

    fHardwareMap[ch] = 1000*mf + 100*mfc + hub;
    printf("Adding ch %i -> %i %i %i\n", ch, mf, mfc, hub);
  }


  std::string line;
  std::getline(f, line);
  std::istringstream ss;

  for ( ; std::getline(f, line); ) {
    ss.clear();
    ss.str(line.c_str());
    ch = row = col = 0;
    ss >> ch >> roc >> col >> row;

    // Just remember that on the plane tester it's channel 22

    if (ch  >  MAXCHNS) { printf("ERROR: over MAXCHNS %i\n", ch); };
    if (row >= MAXROWS) { printf("ERROR: over MAXROWS %i\n", row); };
    if (col >= MAXCOLS) { printf("ERROR: over MAXCOLS %i\n", col); };
    if (roc >= MAXROCS) { printf("ERROR: over MAXROCS %i\n", roc); };
    if (PLTGainCal::DEBUGLEVEL) {
      printf("%i %i %i\n", ch, row, col);
    }

    irow = RowIndex(row);
    icol = ColIndex(col);
    ich  = ChIndex(ch);

    if (irow < 0 || icol < 0 || ich < 0) {
      continue;
    }

    ss >> GC[ich][roc][icol][irow][0]
       >> GC[ich][roc][icol][irow][1]
       >> GC[ich][roc][icol][irow][2]
       >> GC[ich][roc][icol][irow][3]
       >> GC[ich][roc][icol][irow][4];

    // dude, you really don't want to do this..
    if (PLTGainCal::DEBUGLEVEL) {
      for (int i = 0; i != NROCS; ++i) {
        for (int j = 0; j != 5; ++j) {
          printf("%6.2E ", GC[ich][i][icol][irow][j]);
        }
      }
      printf("\n");
    }

  }

  // Apparently this file was read no problem...
  fIsGood = true;


  return;
}


void PLTGainCal::ReadGainCalFileExt (std::string const GainCalFileName, int const roc)
{
  int const ch = 1;
  int row, col;
  int irow;
  int icol;
  int ich;

  //int const mf = 8, mfc = 1, hub = 5;
  fHardwareMap[ch] = 1000*ch + roc;// DA: TODO unused?
  printf("Adding ch %i as -> %i\n", ch, fHardwareMap[ch]);

  ifstream f(GainCalFileName.c_str());
  if (!f) {
    std::cerr << "ERROR: cannot open file: " << GainCalFileName << std::endl;
    throw;
  }

  // Loop over header lines in the input data file
  TString FunctionLine;
  FunctionLine.ReadLine(f);
  FunctionLine.ReadLine(f);
  FunctionLine.ReplaceAll("par[", "[");

  // Set the root function
  TF1 MyFunction("GainCalFitFunction", FunctionLine, -10000, 10000);
  fFitFunction = MyFunction;

  // Get blank line out of the way
  FunctionLine.ReadLine(f);

  std::string line;
  std::getline(f, line);
  std::istringstream ss;

  // extra word in gaincalfile
  std::string PixWord;

  // Just check this number to make sure it's the right size...
  float Coefs[4];
  if (fNParams > 4) {
    std::cerr << "ERROR: NParams is too huge PLTGainCal::ReadGainCalFileExt()" << std::endl;
    exit(1);
  }
  for ( ; std::getline(f, line); ) {
    ss.clear();
    ss.str(line.c_str());

    for (int ipar = 0; ipar < fNParams; ++ipar) {
      ss >> Coefs[ipar];
    }
    ss >> PixWord
       >> col
       >> row;


    if (row >= MAXROWS) { printf("ERROR: over MAXROWS %i\n", row); };
    if (col >= MAXCOLS) { printf("ERROR: over MAXCOLS %i\n", col); };
    if (PLTGainCal::DEBUGLEVEL) {
      printf("%i %i %i\n", ch, row, col);
    }

    irow = RowIndex(row);
    icol = ColIndex(col);
    ich  = ChIndex(ch);

    if (irow < 0 || icol < 0 || ich < 0) {
      continue;
    }

    for (int ipar = 0; ipar < fNParams; ++ipar) {
      GC[ich][roc][icol][irow][ipar] = Coefs[ipar];
    }

    // dude, you really don't want to do this..
    if (PLTGainCal::DEBUGLEVEL) {
      for (int i = 0; i != NROCS; ++i) {
        for (int j = 0; j != 5; ++j) {
          printf("%6.2E ", GC[ich][i][icol][irow][j]);
        }
      }
      printf("\n");
    }

  }

  // Apparently this file was read no problem...
  fIsGood = true;


  return;
}


void PLTGainCal::CheckGainCalFile(std::string const GainCalFileName, int const Channel)
{
  ReadGainCalFile(GainCalFileName, 0);

  int const ich  = ChIndex(Channel);

  int NMissing = 0;
  int NTotal = 0;

  for (int j = 0; j != NROCS; ++j) {
    for (int k = 0; k != PLTU::NCOL; ++k) {
      for (int m = 0; m != PLTU::NROW; ++m) {
        ++NTotal;
        if (
          GC[ich][j][k][m][0] == 0 &&
          GC[ich][j][k][m][1] == 0 &&
          GC[ich][j][k][m][2] == 0 &&
          GC[ich][j][k][m][3] == 0 &&
          GC[ich][j][k][m][4] == 0) {
          printf("Missing Coefs: iCh %2i  iRoc %1i  iCol %2i  iRow %2i\n", ich, j, k, m);
          ++NMissing;
        }
      }
    }
  }

  printf("Number missing is: %i / %i = %E\n", NMissing, NTotal, (float) NMissing / (float) NTotal);

  return;
}


void PLTGainCal::PrintGainCal5 ()
{
  // dude, you really don't want to do this..
  for (int ich = 0; ich != MAXCHNS; ++ich) {
    for (int iroc = 0; iroc != NROCS; ++iroc) {
      for (int icol = 0; icol != 52; ++icol) {
        for (int irow = 0; irow != 80; ++irow) {

          for (int j = 0; j != 5; ++j) {
            printf("%6.2E ", GC[ich][iroc][icol][irow][j]);
          }
          printf("\n");
        }
      }
    }
  }

  return;
}


void PLTGainCal::PrintGainCal (FILE* f)
{
  // dude, you really don't want to do this..
  for (int ich = 0; ich != MAXCHNS; ++ich) {
    for (int iroc = 0; iroc != NROCS; ++iroc) {
      for (int icol = 0; icol != 52; ++icol) {
        for (int irow = 0; irow != 80; ++irow) {

          if (f) {
            fprintf(f, "%2i %2i %2i %2i ", ich+1, iroc, icol, irow);
          } else {
            printf("%2i %2i %2i %2i ", ich+1, iroc, icol, irow);
          }
          for (int j = 0; j != fNParams; ++j) {
            if (f) {
              fprintf(f, "%15.6E ", GC[ich][iroc][icol][irow][j]);
            } else {
              printf("%15.6E ", GC[ich][iroc][icol][irow][j]);
            }
          }
          if (f) {
            fprintf(f, "\n");
          } else {
            printf("\n");
          }
        }
      }
    }
  }

  return;
}

void PLTGainCal::ReadGainCalFile3 (std::string const GainCalFileName)
{
  int mFec, mFecChannel, hubAddress;
  int ch(1), row, col, roc;
  int irow;
  int icol;
  int ich;
  int iroc;

  ifstream f(GainCalFileName.c_str());
  if (!f) {
    std::cerr << "ERROR: cannot open file: " << GainCalFileName << std::endl;
    throw;
  }



  // Loop over header lines in the input data file
  for (std::string line; std::getline(f, line); ) {
    if (line == "") {
      break;
    }
  }

  std::istringstream ss;
  for (std::string line ; std::getline(f, line); ) {
    ss.clear();
    ss.str(line.c_str());
    ss >> mFec >> mFecChannel >> hubAddress >> roc >> col >> row;

    if (ch  >= MAXCHNS) { printf("ERROR: over MAXCHNS\n"); };
    if (row >= MAXROWS) { printf("ERROR: over MAXROWS\n"); };
    if (col >= MAXCOLS) { printf("ERROR: over MAXCOLS\n"); };
    if (PLTGainCal::DEBUGLEVEL) {
      printf("%i %i %i\n", ch, row, col);
    }

    irow = RowIndex(row);
    icol = ColIndex(col);
    ich  = ChIndex(ch);
    iroc = RocIndex(roc);

    if (irow < 0 || icol < 0 || ich < 0) {
      continue;
    }

    ss >> GC[ich][iroc][icol][irow][0]
       >> GC[ich][iroc][icol][irow][1]
       >> GC[ich][iroc][icol][irow][2];

    // dude, you really don't want to do this..
    if (PLTGainCal::DEBUGLEVEL) {
      for (int i = 0; i != NROCS; ++i) {
        printf("%1i %1i %2i %1i %2i %2i", mFec, mFecChannel, hubAddress, roc, col, row);
        for (int j = 0; j != 3; ++j) {
          printf(" %9.1E", GC[ich][i][icol][irow][j]);
        }
        printf("\n");
      }
    }
  }

  // Apparently this file was read no problem...
  fIsGood = true;

  return;
}



void PLTGainCal::ReadTesterGainCalFile (std::string const GainCalFileName)
{
  // Set NParams for 3
  std::cout << "Reading file for Tester GainCal and setting NParams to 3: " << GainCalFileName << std::endl;
  fNParams = 3;

  // Zero for params that don't exist for the teststand
  int const mFec = 0;
  int const mFecChannel = 0;
  int const hubAddress = 0;
  int const ch = 1;
  int const roc = 0;
  int row, col;
  int irow;
  int icol;
  int ich;
  int iroc;


  // If you supply a blank name you did so on purpose (or should have!!)
  if (GainCalFileName == "") {
    std::cout << "ReadTesterGainCalFile sees you have blank input filename.  All gain coefficients are set to zero." << std::endl;
    fIsGood = true;
    return;
  }


  // Open the gaincal file.  If I can't read it we are in trouble.. throw something ugly.
  ifstream f(GainCalFileName.c_str());
  if (!f) {
    std::cerr << "ERROR: cannot open file: " << GainCalFileName << std::endl;
    throw;
  }


  std::istringstream ss;
  for (std::string line ; std::getline(f, line); ) {
    ss.clear();
    ss.str(line.c_str());
    ss >> col >> row;

    if (ch  >= MAXCHNS) { printf("ERROR: over MAXCHNS\n"); };
    if (row >= MAXROWS) { printf("ERROR: over MAXROWS\n"); };
    if (col >= MAXCOLS) { printf("ERROR: over MAXCOLS\n"); };
    if (PLTGainCal::DEBUGLEVEL) {
      printf("%i %i %i\n", ch, row, col);
    }

    irow = RowIndex(row);
    icol = ColIndex(col);
    ich  = ChIndex(ch);
    iroc = RocIndex(roc);

    if (irow < 0 || icol < 0 || ich < 0) {
      continue;
    }

    ss >> GC[ich][iroc][icol][irow][0]
       >> GC[ich][iroc][icol][irow][1]
       >> GC[ich][iroc][icol][irow][2];

    // dude, you really don't want to do this..
    if (PLTGainCal::DEBUGLEVEL) {
      for (int i = 0; i != NROCS; ++i) {
        printf("%1i %1i %2i %1i %2i %2i", mFec, mFecChannel, hubAddress, roc, col, row);
        for (int j = 0; j != 3; ++j) {
          printf(" %9.1E", GC[ich][i][icol][irow][j]);
        }
        printf("\n");
      }
    }
  }

  // Apparently this file was read no problem...
  fIsGood = true;

  return;
}



void PLTGainCal::ResetGC ()
{
  GC.resize(0);

  // Reset everything
  for (int i = 0; i != NCHNS; ++i) {

    std::vector<std::vector<std::vector<std::vector<float> > > > tmp_4d;

    for (int j = 0; j != NROCS; ++j) {

      std::vector<std::vector<std::vector<float> > > tmp_3d;

      for (int k = 0; k != PLTU::NCOL; ++k) {

	std::vector<std::vector<float> > tmp_2d;

        for (int m = 0; m != PLTU::NROW; ++m) {

	  std::vector<float> tmp_1d;
          for (int n = 0; n != 5; ++n) {
	    tmp_1d.push_back(0.);
          }

	tmp_2d.push_back(tmp_1d);
        }
      tmp_3d.push_back(tmp_2d);
      }
    tmp_4d.push_back(tmp_3d);
    }
  GC.push_back(tmp_4d);
  }
  return;
}
