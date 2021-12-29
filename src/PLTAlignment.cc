#include "PLTAlignment.h"

#include <map>
#include <cstdlib>
#include "GetNames.h"
#include "Utils.h"

using namespace std;


PLTAlignment::PLTAlignment (): ErrorsFromFile(false), fIsGood(false) {

  for (int i=0; i != 7; i++){// DA: TODO make error bars equal for all DUTs in X and Y
    fErrorsX.push_back(0.015);
    fErrorsY.push_back(0.015);
  }
}


void PLTAlignment::ReadAlignmentFile(const string & in_file_name, bool use_raw) {

  // So far so good..
  fIsGood = true;
  fTelescopeMap.clear();
  /** if use raw alignmnent, use raw telescope id */
  int16_t telescope_id(tel::Config::telescope_id_);
  if (use_raw) {
    telescope_id = GetRawID();
    tel::warning(Form("Using raw alignment with telescope ID %i", telescope_id));
  }

  // Open file
  std::cout << "Opening " << in_file_name << std::endl;
  std::ifstream InFile(in_file_name.c_str());
  if (!InFile.is_open()) {
    fIsGood = false;
    std::cerr << "ERROR: cannot open alignment constants filename: " << in_file_name << std::endl;
    throw;
  }

  // Read each line in file
  int tid, channel, roc;
  for (std::string InLine; std::getline(InFile, InLine); ) {
    if (InLine.size() < 3 or InLine.at(0) == '#') { continue; }
    std::istringstream LineStream;
    LineStream.str(InLine);

    LineStream >> tid >> channel >> roc;
    if (tid != telescope_id) { continue; }  // only fill the data for the correct telescope id

    std::pair<int, int> CHROC = std::make_pair(channel, roc);

    // read one line at a time
    float R, RZ, RY, X, Y, Z, dX(0), dY(0);

    // If the ROC is -1 it is telescope coords, 0,1,2 are ROCs, anything else is bad.
    if (roc == 0 and fTelescopeMap.count(channel) == 0) { // use 0 if roc -1 is not in the data
      fTelescopeMap[channel].GRZ = 0;
      fTelescopeMap[channel].GRY = 0;
      fTelescopeMap[channel].GX  = 0;
      fTelescopeMap[channel].GY  = 0;
      fTelescopeMap[channel].GZ  = 0;
    }
    if (roc == -1) {
      LineStream >> RZ >> RY >> X >> Y >> Z;
      fTelescopeMap[channel].GRZ = RZ;
      fTelescopeMap[channel].GRY = RY;
      fTelescopeMap[channel].GX  = X;
      fTelescopeMap[channel].GY  = Y;
      fTelescopeMap[channel].GZ  = Z;
    } else if (roc < 7) {

      LineStream >> R >> X  >> Y >> Z >> dX >> dY;

      // Construct the alignment obj
      CP C;
      C.LR = R;
      C.LX = X;
      C.LY = Y;
      C.LZ = Z;
      C.GRZ = fTelescopeMap[channel].GRZ;
      C.GRY = fTelescopeMap[channel].GRY;
      C.GX = fTelescopeMap[channel].GX;
      C.GY = fTelescopeMap[channel].GY;
      C.GZ = fTelescopeMap[channel].GZ;

      // Save this to alignment map
      fConstantMap[CHROC] = C;
    } else {
      std::cerr << "WARNING: Alignment file contains things I do not recognize: " << InLine << std::endl;
    }

    /** try to read the errors from file */
    if (dX > 0 and dY > 0){
      ErrorsFromFile = true;
      SetErrorX(roc, dX);
      SetErrorY(roc, dY);
    }

  }

  // Close input file
  InFile.close();

  if (fTelescopeMap.empty()){
    tel::warning(Form("Did not find telescope %i in the alignments file %s", tel::Config::telescope_id_, tel::split(in_file_name, '/').back().c_str()));
    ReadAlignmentFile(in_file_name, true);
  }
} // end ReadAlignmentFile


string PLTAlignment::GetAlignment(uint16_t telescope_id, uint16_t n_rocs, bool write_errors) {

  ostringstream os;
  for (auto it: fTelescopeMap){
    int const channel = it.first;
    TelescopeAlignmentStruct & tel = it.second;
    if (tel.GRZ != 0){
      os << Form("% 3i% 4i  -1 %+1.4E  %+1.4E  %+1.4E  %+1.4E  %+1.4E", telescope_id, channel, tel.GRZ, tel.GRY, tel.GX, tel.GY, tel.GZ) << endl;
    }
    for (int iroc = 0; iroc < n_rocs; iroc++){
      pair<int, int> ch_roc = make_pair(channel, iroc);
      if (fConstantMap.count(ch_roc) == 0) {
        std::cerr << "ERROR: No entry in fConstantMap for Ch ROC: " << channel << " " << iroc << endl;
        continue;
      }
      CP & C = fConstantMap[ch_roc];
      os << Form("% 3i% 4i% 4i  %+1.4E  %+24.4E  %+1.4E  %+1.4E", telescope_id, channel, iroc, C.LR, C.LX, C.LY, C.LZ);
      if (write_errors) { os << Form("  %+1.4E  %+1.4E", GetErrorX(iroc), GetErrorY(iroc)); }
      os << endl;
    }
  }
  return os.str();
}


void PLTAlignment::WriteAlignmentFile (const uint16_t telescop_id, const uint16_t n_rocs, bool write_errors) {
  /** read alignmnent file, overwrite settings if already existing, create new entry if not existing */
  // get lines from file
  vector<string> lines;
  string line_str;
  ifstream in(GetAlignmentFilename());
  if (!in) {
    std::cerr << "ERROR: cannot open file: " << GetAlignmentFilename() << endl;
    throw;
  }
  cout << "Writing aligment to file " << GetAlignmentFilename() << endl;
  while (getline(in, line_str)) { lines.emplace_back(line_str); }
  in.close();
  // write new settings
  ofstream out(GetAlignmentFilename());
  bool wrote_data(false);
  for (const auto & line: lines) {
    if (line.find('#') != string::npos) {
      out << line << endl;
      continue; }
    istringstream s(line);
    int tel, ch, roc;
    float rz, x, y, z, dx(0), dy(0);
    s >> tel >> ch >> roc >> rz >> x >> y >> z >> dx >> dy;
    if (tel != telescop_id) { out << line << endl; }
    else {
      SetErrorX(roc, dx); SetErrorY(roc, dy);
      if (not wrote_data and roc == n_rocs - 1){
        out << GetAlignment(telescop_id, n_rocs, (dx > 0 and dy > 0) or write_errors); // we have an existing alignment
        wrote_data = true;
      }
    }
  }
  if (not wrote_data) {  // we have a new alignment
    const int comment_length = 88;
    out << Form("# TELESCOPE % 2i ", tel::Config::telescope_id_) << string(comment_length, '-') << endl;
    out << GetAlignment(telescop_id, n_rocs, write_errors);
  }
  out.close();
}


float PLTAlignment::PXtoLX (int const px)
{
  return PLTU::PIXELWIDTH * ((px + 0.0000000001) - PLTU::DIACENTERX);
}

float PLTAlignment::PYtoLY (int const py)
{
  return PLTU::PIXELHEIGHT * ((py + 0.0000000001) - PLTU::DIACENTERY);
}

std::pair<int, int> PLTAlignment::PXYfromLXY (std::pair<float, float> const& LXY)
{
  return std::make_pair( PXfromLX(LXY.first), PYfromLY(LXY.second));
}

std::pair<float, float> PLTAlignment::PXYDistFromLXYDist (std::pair<float, float> const& LXYDist)
{
  return std::make_pair( LXYDist.first / (float)  PLTU::PIXELWIDTH, LXYDist.second / (float)  PLTU::PIXELHEIGHT);
}


void PLTAlignment::AlignHit (PLTHit& Hit)
{
  // Grab the constants and check that they are there..
  CP* C = GetCP(Hit.Channel(), Hit.ROC());
  if (C == 0x0) {
    std::cerr << "ERROR: This is not in the aligment constants map: Channel:" << Hit.Channel() << "  ROC:" << Hit.ROC() << std::endl;
    return;
  }

  int const PX = Hit.Column();
  int const PY = Hit.Row();

  // set w.r.t. center of diamond
  float LX = PXtoLX(PX);
  float LY = PYtoLY(PY);

//  std::cout << LX << " " << LY << std::endl;

  std::vector<float> TXYZ;
  LtoTXYZ(TXYZ, LX, LY, Hit.Channel(), Hit.ROC());


  if (DEBUG) {
    printf("TtoL - L XY DIFF %12.3f %12.3f\n",
        TtoLX(TXYZ[0], TXYZ[1], Hit.Channel(), Hit.ROC()) - LX,
        TtoLY(TXYZ[0], TXYZ[1], Hit.Channel(), Hit.ROC()) - LY);
  }

  std::vector<float> GXYZ;
  TtoGXYZ(GXYZ, TXYZ[0], TXYZ[1], TXYZ[2], Hit.Channel(), Hit.ROC());


  // Set the local, telescope, and global hit coords
  Hit.SetLXY(LX, LY);
  Hit.SetTXYZ(TXYZ[0], TXYZ[1], TXYZ[2]);
  Hit.SetGXYZ(GXYZ[0], GXYZ[1], GXYZ[2]);

  //printf("Channel %2i ROC %1i  Col %2i Row %2i  %12.3E  %12.3E - %12.3E  %12.3E  %12.3E - %12.3E  %12.3E  %12.3E\n",
  //    Hit.Channel(), Hit.ROC(), Hit.Column(), Hit.Row(), LX, LY, TXYZ[0], TXYZ[1], TXYZ[2], GXYZ[0], GXYZ[1], GXYZ[2]);

  //std::vector<float> TV;
  //GtoTXYZ(TV, GX, GY, GZ, Hit.Channel(), Hit.ROC());
  //if (DEBUG) {
  //  printf("GtoT - T XYZ DIFF %12.3f %12.3f %12.3f\n",
  //      TV[0] - TX, TV[1] - TY, TV[2] - TZ);
  //}

  return;
}



float PLTAlignment::GetTZ (int const Channel, int const ROC)
{
  return GetCP(Channel, ROC)->LZ;
}




float PLTAlignment::TtoLX (float const TX, float const TY, int const Channel, int const ROC)
{
  return TtoLXY(TX, TY, Channel, ROC).first;
}


float PLTAlignment::TtoLY (float const TX, float const TY, int const Channel, int const ROC)
{
  return TtoLXY(TX, TY, Channel, ROC).second;
}


std::pair<float, float> PLTAlignment::TtoLXY (float const TX, float const TY, int const Channel, int const ROC)
{
  std::pair<int, int> CHROC = std::make_pair(Channel, ROC);
  CP* C = fConstantMap.count(CHROC) == 1 ? &fConstantMap[CHROC] : (CP*) 0x0;

  if (!C) {
//    std::cerr << "ERROR: cannot grab the constant mape for this CH ROC: " << CHROC.first << " " << CHROC.second << std::endl;
    return std::make_pair(-999, -999);
  }


  float const LXA = TX - C->LX;
  float const LYA = TY - C->LY;


  float const cl = cos(C->LR);
  float const sl = sin(C->LR);


  float const LX =  LXA * cl + LYA * sl;
  float const LY = -LXA * sl + LYA * cl;

  //printf("XY DIFF %12.3f  %12.3f\n", TX - LX, TY - LY);

  return std::make_pair(LX, LY);
}


void PLTAlignment::VTtoVGXYZ (std::vector<float>& VOUT, float const TX, float const TY, float const TZ, int const Channel, int const ROC)
{
  // Get the constants for this telescope/plane etc
  std::pair<int, int> CHROC = std::make_pair(Channel, ROC);
  CP* C = fConstantMap.count(CHROC) == 1 ? &fConstantMap[CHROC] : (CP*) 0x0;

  if (!C) {
    std::cerr << "ERROR: cannot grab the constant mape for this CH ROC: " << CHROC.first << " " << CHROC.second << std::endl;
    return;
    throw;
  }

  // Global rotation about Zaxis
  float const cgz = cos(C->GRZ);
  float const sgz = sin(C->GRZ);
  float GXZ = TX * cgz - TY * sgz;
  float GYZ = TX * sgz + TY * cgz;
  float GZZ = TZ;


  // Global rotation about Yaxis
  float const cgy = cos(C->GRY);
  float const sgy = sin(C->GRY);
  float GXY = GXZ * cgy + GZZ * sgy;
  float GYY = GYZ;
  float GZY = GZZ * cgy - GXZ * sgy;

  VOUT.resize(3, 0);
  VOUT[0] = GXY;
  VOUT[1] = GYY;
  VOUT[2] = GZY;


  return;
}

void PLTAlignment::LtoTXYZ (std::vector<float>& VOUT, float const LX, float const LY, int const Channel, int const ROC)
{
  std::pair<int, int> CHROC = std::make_pair(Channel, ROC);
  CP* C = fConstantMap.count(CHROC) == 1 ? &fConstantMap[CHROC] : (CP*) 0x0;
//  std::cout << C->LR << std::endl;

  if (!C) {
    std::cerr << "ERROR: cannot grab the constant map for this CH ROC: " << CHROC.first << " " << CHROC.second << std::endl;
    throw;
  }

  // Start with local rotation
  float const cl = cos(C->LR);
  float const sl = sin(C->LR);
  float TX = LX * cl - LY * sl;
  float TY = LX * sl + LY * cl;
  float TZ = C->LZ;

  // Local translation
  TX += C->LX;
  TY += C->LY;

  VOUT.resize(3);
  VOUT[0] = TX;
  VOUT[1] = TY;
  VOUT[2] = TZ;

  return;
}


void PLTAlignment::TtoGXYZ (std::vector<float>& VOUT, float const TX, float const TY, float const TZ, int const Channel, int const ROC)
{
  // Get the constants for this telescope/plane etc
  std::pair<int, int> CHROC = std::make_pair(Channel, ROC);
  CP* C = fConstantMap.count(CHROC) == 1 ? &fConstantMap[CHROC] : (CP*) 0x0;

  if (!C) {
    std::cerr << "ERROR: cannot grab the constant mape for this CH ROC: " << CHROC.first << " " << CHROC.second << std::endl;
    return;
    throw;
  }

  // Global rotation about Zaxis
  float const cgz = cos(C->GRZ);
  float const sgz = sin(C->GRZ);
  float GXZ = TX * cgz - TY * sgz;
  float GYZ = TX * sgz + TY * cgz;
  float GZZ = TZ;

  GXZ += C->GX;
  GYZ += C->GY;
  GZZ += C->GZ;


  // Global rotation about Yaxis
  float const cgy = cos(C->GRY);
  float const sgy = sin(C->GRY);
  float GXY = GXZ * cgy + GZZ * sgy;
  float GYY = GYZ;
  float GZY = GZZ * cgy - GXZ * sgy;

  VOUT.resize(3, 0);
  VOUT[0] = GXY;
  VOUT[1] = GYY;
  VOUT[2] = GZY;


  return;
}



void PLTAlignment::LtoGXYZ (std::vector<float>& VOUT, float const LX, float const LY, int const Channel, int const ROC)
{
  std::vector<float> T;
  LtoTXYZ(T, LX, LY, Channel, ROC);
  TtoGXYZ(VOUT, T[0], T[1], T[2], Channel, ROC);
}



void PLTAlignment::GtoTXYZ (std::vector<float>& VOUT, float const GX, float const GY, float const GZ, int const Channel, int const ROC)
{
  // This translates global coordinates back to the telescope coordinates

  // Get the constants for this telescope/plane etc
  std::pair<int, int> CHROC = std::make_pair(Channel, ROC);
  CP* C = fConstantMap.count(CHROC) == 1 ? &fConstantMap[CHROC] : (CP*) 0x0;

  if (!C) {
    std::cerr << "ERROR: cannot grab the constant mape for this CH ROC: " << CHROC.first << " " << CHROC.second << std::endl;
    return;
    throw;
  }

  // Global rotation about Yaxis
  float const cgy = cos(C->GRY);
  float const sgy = sin(C->GRY);
  float GXY = GX * cgy - GZ * sgy;
  float GYY = GY;
  float GZY = GZ * cgy + GX * sgy;


  float const GXA = GXY - C->GX;
  float const GYA = GYY - C->GY;
  float const GZA = GZY - C->GZ;


  float const cl = cos(C->GRZ);
  float const sl = sin(C->GRZ);


  float const TX =  GXA * cl + GYA * sl;
  float const TY = -GXA * sl + GYA * cl;
  float const TZ =  GZA;


  VOUT.resize(3);
  VOUT[0] = TX;
  VOUT[1] = TY;
  VOUT[2] = TZ;

  return;
}








float PLTAlignment::LR (int const ch, int const roc)
{
  return fConstantMap[ std::make_pair(ch, roc) ].LR;
}


float PLTAlignment::LX (int const ch, int const roc)
{
  return fConstantMap[ std::make_pair(ch, roc) ].LX;
}


float PLTAlignment::LY (int const ch, int const roc)
{
  return fConstantMap[ std::make_pair(ch, roc) ].LY;
}


float PLTAlignment::LZ (int const ch, int const roc)
{
  return fConstantMap[ std::make_pair(ch, roc) ].LZ;
}


float PLTAlignment::GRZ (int const ch, int const roc)
{
  return fConstantMap[ std::make_pair(ch, roc) ].GRZ;
}

float PLTAlignment::GRY (int const ch, int const roc)
{
  return fConstantMap[ std::make_pair(ch, roc) ].GRY;
}


float PLTAlignment::GX (int const ch, int const roc)
{
  return fConstantMap[ std::make_pair(ch, roc) ].GX;
}


float PLTAlignment::GY (int const ch, int const roc)
{
  return fConstantMap[ std::make_pair(ch, roc) ].GY;
}


float PLTAlignment::GZ (int const ch, int const roc)
{
  return fConstantMap[ std::make_pair(ch, roc) ].GZ;
}

PLTAlignment::CP* PLTAlignment::GetCP (int const ch, int const roc)
{
  std::pair<int, int> CHROC = std::make_pair(ch, roc);
  if (fConstantMap.count(CHROC)) {
    return &fConstantMap[ CHROC ];
  }

  return (CP*) 0x0;
}

PLTAlignment::CP* PLTAlignment::GetCP (std::pair<int, int> const& CHROC)
{
  if (fConstantMap.count(CHROC)) {
    return &fConstantMap[ CHROC ];
  }

  return (CP*) 0x0;
}

std::vector< std::pair<int, int> > PLTAlignment::GetListOfChannelROCs ()
{
  // returns a vector containing the pixel fed channels for everything in the alignment

  std::vector< std::pair<int, int> > ROCS;
  for (std::map< std::pair<int, int>, CP >::iterator i = fConstantMap.begin(); i != fConstantMap.end(); ++i) {
    ROCS.push_back(i->first);
  }

  return ROCS;
}

std::vector<int> PLTAlignment::GetListOfChannels ()
{
  std::vector<int> Channels;

  for (std::map<int, TelescopeAlignmentStruct>::iterator it = fTelescopeMap.begin(); it != fTelescopeMap.end(); ++it)
  {
    Channels.push_back(it->first);
  }

  return Channels;
}

void PLTAlignment::AddToLR (int const ch, int const roc, float val)
{
  float oldval = fConstantMap[ std::make_pair(ch, roc) ].LR;
  fConstantMap[ std::make_pair(ch, roc) ].LR = oldval+val;
}


void PLTAlignment::AddToLX (int const ch, int const roc, float val)
{
  float oldval = fConstantMap[ std::make_pair(ch, roc) ].LX;
  fConstantMap[ std::make_pair(ch, roc) ].LX = oldval+val;
}


void PLTAlignment::AddToLY (int const ch, int const roc, float val)
{
  float oldval = fConstantMap[ std::make_pair(ch, roc) ].LY;
  fConstantMap[ std::make_pair(ch, roc) ].LY = oldval+val;
}


void PLTAlignment::AddToLZ (int const ch, int const roc, float val)
{
  float oldval = fConstantMap[ std::make_pair(ch, roc) ].LZ;
  fConstantMap[ std::make_pair(ch, roc) ].LZ = oldval+val;
}

void PLTAlignment::AddToGX (int const ch, float val)
{
  float oldval = fTelescopeMap[ch].GX;
  fTelescopeMap[ch].GX  = oldval+val;
  fConstantMap[ std::make_pair(ch, 0) ].GX = oldval+val;
  fConstantMap[ std::make_pair(ch, 1) ].GX = oldval+val;
  fConstantMap[ std::make_pair(ch, 2) ].GX = oldval+val;
}


void PLTAlignment::AddToGY (int const ch, float val)
{
  float oldval = fTelescopeMap[ch].GY;
  fTelescopeMap[ch].GY  = val;
  fConstantMap[ std::make_pair(ch, 0) ].GY = oldval+val;
  fConstantMap[ std::make_pair(ch, 1) ].GY = oldval+val;
  fConstantMap[ std::make_pair(ch, 2) ].GY = oldval+val;
}


void PLTAlignment::AddToGZ (int const ch, float val)
{
  float oldval = fTelescopeMap[ch].GZ;
  fTelescopeMap[ch].GZ  = val;
  fConstantMap[ std::make_pair(ch, 0) ].GZ = oldval+val;
  fConstantMap[ std::make_pair(ch, 1) ].GZ = oldval+val;
  fConstantMap[ std::make_pair(ch, 2) ].GZ = oldval+val;
}

void PLTAlignment::SetErrors(int telescopeID, bool initial){

    uint8_t id = telescopeID;
    uint8_t n_planes = GetNPlanes();
    if (initial) { /** just use increased digital resolution for initial config */
      for (auto i_plane(0); i_plane < n_planes; i_plane++){
        SetErrorX(i_plane, PLTU::PIXELWIDTH / sqrt(12) * 2.5);
        SetErrorY(i_plane, PLTU::PIXELHEIGHT / sqrt(12) * 2.5);
      }
    } else if (ErrorsFromFile) {
      std::cout << "Already found errors in alignment file!" << std::endl;
    } else {
        if (telescopeID == 1){
            SetErrorX( 0, 0.0138452);
            SetErrorX( 1, 0.0041562);
            SetErrorX( 2, 0.00850018);
            SetErrorX( 3, 0.00796671);
            SetErrorX( 4, 0.00778308);
            SetErrorX( 5, 0.0142984);

            SetErrorY( 0, 0.013026);
            SetErrorY( 1, 0.00384812);
            SetErrorY( 2, 0.00728567);
            SetErrorY( 3, 0.00688992);
            SetErrorY( 4, 0.00460671);
            SetErrorY( 5, 0.0130273);
        }
        else if (telescopeID == 2){
            SetErrorX( 0, 0.0116691);
            SetErrorX( 1, 0.00476967);
            SetErrorX( 2, 0.00974399);
            SetErrorX( 3, 0.00772797);
            SetErrorX( 4, 0.00423922);
            SetErrorX( 5, 0.0116691);

            SetErrorY( 0, 0.0135081);
            SetErrorY( 1, 0.00424304);
            SetErrorY( 2, 0.00661225);
            SetErrorY( 3, 0.00689123);
            SetErrorY( 4, 0.00310896);
            SetErrorY( 5, 0.0104977);
        }
        else if (telescopeID == 5){
            SetErrorX( 0, 0.01);
            SetErrorX( 1, 0.01);
            SetErrorX( 2, 0.01);
            SetErrorX( 3, 0.01);

            SetErrorY( 0, 0.01);
            SetErrorY( 1, 0.01);
            SetErrorY( 2, 0.01);
            SetErrorY( 3, 0.01);

        }
        else if (telescopeID == 7){
            SetErrorX( 0, 0.01);
            SetErrorX( 1, 0.01);
            SetErrorX( 2, 0.01);
            SetErrorX( 3, 0.01);

            SetErrorY( 0, 0.01);
            SetErrorY( 1, 0.01);
            SetErrorY( 2, 0.01);
            SetErrorY( 3, 0.01);
        }
        else if (id == 9 || id == 11 || id == 12 || id>=14){
            float correction = 2.5;
            bool use_geom_fac = false;
            float geom_fac = use_geom_fac ? sqrt(5) : 1;
            for (uint8_t i = 0; i < 4; i++){
                if (i == 0 or i == 3) {
                    SetErrorX(i, geom_fac * 0.015 / sqrt(12) * correction);
                    SetErrorY(i, geom_fac * 0.01 / sqrt(12) * correction);
                }
                else {
                    SetErrorX(i, 0.015/sqrt(12) * correction);
                    SetErrorY(i, 0.01/sqrt(12) * correction);
                }
            }
//            SetErrorX( 0, sqrt(5)*0.015/sqrt(12));
//            SetErrorX( 1, 0.015/sqrt(12));
//            SetErrorX( 2, 0.015/sqrt(12));
//            SetErrorX( 3, sqrt(5)*0.015/sqrt(12));
//
//            SetErrorY( 0, sqrt(5)*0.01/sqrt(12));
//            SetErrorY( 1, 0.01/sqrt(12));
//            SetErrorY( 2, 0.01/sqrt(12));
//            SetErrorY( 3, sqrt(5)*0.01/sqrt(12));
        }
        else if (telescopeID == 8){
            SetErrorX( 0, 0.01);
            SetErrorX( 1, 0.01);
            SetErrorX( 2, 0.01);
            SetErrorX( 3, 0.01);
            SetErrorX( 4, 0.01);
            SetErrorX( 5, 0.01);

            SetErrorY( 0, 0.01);
            SetErrorY( 1, 0.01);
            SetErrorY( 2, 0.01);
            SetErrorY( 3, 0.01);
            SetErrorY( 4, 0.01);
            SetErrorY( 5, 0.01);
        }
        else if (id == 10 || id == 13 || id==15){
            SetErrorX( 0, 0.015);
            SetErrorX( 1, 0.015);
            SetErrorX( 2, 0.015);
            SetErrorX( 3, 0.015);
            SetErrorX( 4, 0.015);
            SetErrorX( 5, 0.015);
            SetErrorX( 6, 0.015);

            SetErrorY( 0, 0.01);
            SetErrorY( 1, 0.01);
            SetErrorY( 2, 0.01);
            SetErrorY( 3, 0.01);
            SetErrorY( 4, 0.01);
            SetErrorY( 5, 0.01);
            SetErrorY( 6, 0.01);
        }
        else if (telescopeID == -1){
            SetErrorX( 0, 0.01);
            SetErrorX( 1, 0.01);
            SetErrorX( 2, 0.01);
            SetErrorX( 3, 0.01);
            SetErrorX( 4, 0.01);
            SetErrorX( 5, 0.01);

            SetErrorY( 0, 0.01);
            SetErrorY( 1, 0.01);
            SetErrorY( 2, 0.01);
            SetErrorY( 3, 0.01);
            SetErrorY( 4, 0.01);
            SetErrorY( 5, 0.01);
        }
        else{
            std::cout << "ERROR: No Errors defined for telescopeID==" << telescopeID << std::endl;
            std::cout << "Exiting.." << std::endl;
            std::exit(0);
        }
    }

}

void PLTAlignment::ResetPlane(int const ch, int const roc) {

  fConstantMap[std::make_pair(ch, roc)].LX = 0;
  fConstantMap[std::make_pair(ch, roc)].LY = 0;
  fConstantMap[std::make_pair(ch, roc)].LR = 0;
}
