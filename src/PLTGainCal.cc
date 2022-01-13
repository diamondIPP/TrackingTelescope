#include "PLTGainCal.h"
#include "Utils.h"
#include <algorithm>

#define DEF_CHARGE -9999
#define MAX_VCAL 7 * 255

using namespace std;


PLTGainCal::PLTGainCal () {
  ResetGC();
}

PLTGainCal::PLTGainCal (int nrocs, bool isExternalFunction ): NROCS(nrocs) {
  ResetGC();
  fIsExternalFunction = isExternalFunction;
  ReadVcalCal();
}

PLTGainCal::PLTGainCal (std::string const & GainCalFileName, int const NParams): fNParams(NParams) {
  ResetGC();
  ReadGainCalFile(GainCalFileName);
}


float PLTGainCal::GetCharge(int const ch, int const roc, int const col, int const row, int adc) {
  /** Get charge, note roc number is 0, 1, 2... */
  if (ChIndex(ch)   >= MAXCHNS) { printf("ERROR: over MAXCHNS: %i\n", ch); };
  if (RowIndex(row) >= PLTU::NROW) { printf("ERROR: over MAXROWS: %i\n", row); };
  if (ColIndex(col) >= PLTU::NCOL) { printf("ERROR: over MAXCOLS: %i\n", col); };

  int16_t irow = RowIndex(row), icol = ColIndex(col), ich  = ChIndex(ch), iroc = RocIndex(roc);
  if (irow < 0 || icol < 0 || ich < 0 || iroc < 0) { return DEF_CHARGE; }

  double vcal;

  if (fIsExternalFunction) {  /** external calibration */
    for (int ipar = 0; ipar < fNParams; ++ipar) { fFitFunction.SetParameter(ipar, GC[ich][iroc][icol][irow][ipar]);}
    if (adc + 1 > fFitFunction.GetMaximum() or adc - 1 < fFitFunction.GetMinimum() or adc == 0 and tel::Config::telescope_id_ == 22) { return DEF_CHARGE; }
    vcal = min(max(fFitFunction.GetX(adc), 0.), double(MAX_VCAL));  // contain vcal in range [0, MAX_VCAL]
  }
  else {  /** old calibration */
    if (fNParams == 3) { vcal = float(adc * adc) * GC[ich][iroc][icol][irow][2] + float(adc) * GC[ich][iroc][icol][irow][1] + GC[ich][iroc][icol][irow][0]; }
    else if (fNParams == 5) { vcal = (TMath::Power( (float) adc, 2) * GC[ich][iroc][icol][irow][0] + (float) adc * GC[ich][iroc][icol][irow][1] + GC[ich][iroc][icol][irow][2]
                                   + (GC[ich][iroc][icol][irow][4] != 0 ? TMath::Exp( (adc - GC[ich][iroc][icol][irow][3]) / GC[ich][iroc][icol][irow][4] ) : 0) ); }
    else { tel::critical(Form("ERROR: PLTGainCal::GetCharge() I do not know of that number of fit parameters: %d", fNParams)); exit(1);}
  }

  if (PLTGainCal::DEBUGLEVEL) { printf("%2i %1i %2i %2i %4i %6.1f %8.1f\n", ch, roc, col, row, adc, vcal, VC.at(iroc).first * vcal + VC.at(iroc).second); }

  return VC.at(iroc).first * vcal + VC.at(iroc).second;
}

void PLTGainCal::ReadGainCalFile (const string & GainCalFileName, int roc) {

  if (GainCalFileName.empty()) {
    fIsGood = false;
    return;
  }

  ifstream f(GainCalFileName.c_str());
  if (!f.is_open()) {
    tel::critical(Form("Cannot open gaincal file: %s", GainCalFileName.c_str()));
    throw;
  }

  std::string line;
  getline(f, line);  // read first line
  fIsExternalFunction = line.find("Parameters of the vcal vs. pulse height fits") != string::npos;

  for (; std::getline(f, line); ) { if (line.empty() ) { break; } }  // Loop over header lines in the input data file

  getline(f, line);
  istringstream linestream(line);
  int i = fIsExternalFunction ? 0 : -4;
  for (float junk; linestream >> junk; ++i) { }
  f.close();
  fNParams = i;

  if (fIsExternalFunction) { ReadGainCalFileExt(GainCalFileName, roc); }
  else {
    if (fNParams == 5) { ReadGainCalFile5(GainCalFileName); }
    else if (fNParams == 3) { ReadGainCalFile3(GainCalFileName); }
    else {
      tel::critical("ERROR: I have no idea how many params you have");
      throw;
    }
  }
}

int PLTGainCal::GetHardwareID (int const Channel)
{
  return fHardwareMap[Channel];
}

void PLTGainCal::ReadGainCalFile5 (const string & GainCalFileName)
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


void PLTGainCal::ReadGainCalFileExt (const string & GainCalFileName, int const roc)
{
  int const ch = 1;
  int row, col;
  int irow;
  int icol;
  int ich;

  //int const mf = 8, mfc = 1, hub = 5;
  fHardwareMap[ch] = 1000*ch + roc;// DA: TODO unused?

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
  TF1 MyFunction("GainCalFitFunction", FunctionLine, -MAX_VCAL, MAX_VCAL);
  fFitFunction = MyFunction;
  fFitFunction.SetNpx(180);

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

void PLTGainCal::ReadGainCalFile3 (const string & GainCalFileName)
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

void PLTGainCal::ReadVcalCal() {
  ifstream f(GetDir() + "data/vcal_calibrations.txt");
  int t_id, roc;
  float gain, offset;
  map<pair<int, int>, pair<float, float>> tmp;
  for (string line; getline(f, line);) {
    istringstream s(line);
    s >> t_id >> roc >> gain >> offset;
    tmp[make_pair(t_id, roc)] = make_pair(gain, offset);
  }
  for (uint8_t i_roc(0); i_roc < GetNPlanes(); i_roc++) {
    pair<int, int> id = make_pair(tel::Config::telescope_id_, i_roc);
    pair<int, int> default_id = make_pair(0, i_roc >= tel::Config::n_tel_planes_ and UseDigitalCalibration() ? 0 : 1);
    VC.emplace_back(tmp.at(tmp.find(id) != tmp.end() ? id : default_id));
  }
}
