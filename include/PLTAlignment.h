#ifndef GUARD_PLTAlignment_h
#define GUARD_PLTAlignment_h

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
#include <cmath>

#include "PLTHit.h"
#include "PLTCluster.h"
#include "PLTU.h"

class PLTAlignment
{
  public:
    PLTAlignment ();
    ~PLTAlignment () = default;

    void ReadAlignmentFile(std::string const &in_file_name, bool=false);
    void WriteAlignmentFile (uint16_t, uint16_t, bool=false);
    std::string GetAlignment(uint16_t, uint16_t, bool=false);
    void AlignHit (PLTHit&);
    bool ErrorsFromFile;
    bool IsGood () { return fIsGood; }

    float TtoLX (float , float, int, int);
    float TtoLY (float, float, int, int);
    std::pair<float, float> TtoLXY (float, float, int, int);
    static float LX2PX(float lx) { return lx / PLTU::PIXELWIDTH + PLTU::DIACENTERX; }
    static float LY2PY(float ly) { return ly / PLTU::PIXELHEIGHT + PLTU::DIACENTERY; }


    void LtoTXYZ (std::vector<float>&, float, float, int, int);
    void LtoGXYZ (std::vector<float>&, float, float, int, int);
    void TtoGXYZ (std::vector<float>&, float, float, float, int, int);
    void GtoTXYZ (std::vector<float>&, float, float, float, int, int);
    void VTtoVGXYZ (std::vector<float>&, float, float, float, int, int);
    float GetTZ (int, int);

    float PXtoLX (int);
    float PYtoLY (int);


    static int PXfromLX(float lx) { return int(lx / PLTU::PIXELWIDTH + PLTU::DIACENTERX); }
    static int PYfromLY(float ly) { return int(ly / PLTU::PIXELHEIGHT + PLTU::DIACENTERY); }
    std::pair<int, int> PXYfromLXY (std::pair<float, float> const&);
    std::pair<float, float> PXYDistFromLXYDist (std::pair<float, float> const&);

    // GtoT

    // Need a function for G->T, G->L
    // Need to make T in cm?  G in cm?
    // Function: IsTTrackFiducial(plane)



    //GetPXPY (LX, LY, ch, ROC);

    float LR  (int, int);
    float LX  (int, int);
    float LY  (int, int);
    float LZ  (int, int);
    float GRZ (int, int);
    float GRY (int, int);
    float GX  (int, int);
    float GY  (int, int);
    float GZ  (int, int);

    struct CP {
      // All constants refer to planes as though we were looking FROM the IP
      float LR;  // Local clockwise rotation
      float LX;  // Local X translation
      float LY;  // Local Y translation
      float LZ;  // Local Z translation
      float GRZ; // Global clockwise rotation about Z
      float GRY; // Global clockwise rotation about Y
      float GX;  // Global X translation
      float GY;  // Global Y translation
      float GZ;  // Global Z translation
    };

    CP* GetCP (int, int);
    CP* GetCP (std::pair<int, int> const&);

    std::vector< std::pair<int, int> > GetListOfChannelROCs ();
    std::vector<int> GetListOfChannels ();

    // Mini struct to be used only in reading alignment file
    struct TelescopeAlignmentStruct {
      float GRZ, GRY, GX, GY, GZ;
    };

    void AddToLR (int, int, float);
    void AddToLX (int, int, float);
    void AddToLY (int, int, float);
    void AddToLZ (int, int, float);
    void AddToGX (int, float);
    void AddToGY (int, float);
    void AddToGZ (int, float);

    void ResetPlane(int, int);

    float GetErrorX(int plane){ return fErrorsX[plane];};
    float GetErrorY(int plane){ return fErrorsY[plane];};

    void SetErrorX(int plane, float val ){ fErrorsX[plane]=val;};
    void SetErrorY(int plane, float val ){ fErrorsY[plane]=val;};

    void SetErrors(int telescopeID, bool initial = false);

  private:
    std::map< std::pair<int, int>, CP > fConstantMap;
    std::map<int, TelescopeAlignmentStruct> fTelescopeMap;

    std::vector< float > fErrorsX;
    std::vector< float > fErrorsY;

    bool fIsGood;

    static bool const DEBUG = false;

};









#endif
