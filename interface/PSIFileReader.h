#ifndef GUARD_PSIFileReader_h
#define GUARD_PSIFileReader_h

#include <fstream>
#include <set>

#include "PLTTelescope.h"
#include "PLTGainCal.h"
#include "PSIGainInterpolator.h"
#include "PLTAlignment.h"
#include "PLTTracking.h"

class PSIFileReader : public PLTTelescope, public PLTTracking
{

public:
    PSIFileReader (std::string const, std::string const, int const, bool const, bool const, bool);
    virtual ~PSIFileReader (){};

    virtual bool OpenFile () = 0;
    virtual void ResetFile () = 0 ;

    void Clear ();

    void ReadPixelMask (std::string const);
    void AddToPixelMask( int, int, int, int);
    bool IsPixelMasked (int const);

    void DrawTracksAndHits (std::string const);

    virtual int GetNextEvent () = 0;

    size_t NHits ();
    PLTHit* Hit (size_t);

    PLTGainCal * GetGainCal () { return &fGainCal; }
    PLTAlignment * GetAlignment() { return &fAlignment; }
    const std::set<int> * GetPixelMask(){ return &fPixelMask; }

    long long GetTime () {return fTime;}

    const int  NMAXROCS;


  protected:

    long long fTime;

    std::set<int> fPixelMask;
    std::vector<PLTHit*> fHits;

    PSIGainInterpolator fGainInterpolator;
    PLTGainCal fGainCal;
    PLTAlignment fAlignment;

    std::string fBinaryFileName;
    bool trackOnlyTelescope;

    // Should we use the GainInterpolator instead of GainCal
    // -> Only for Telescope 2 from May 2014 testmeab for now
    bool fUseGainInterpolator;
//    bool trackOnlyTelescope;

    std::map<int, PLTPlane> fPlaneMap;

    std::string fBaseCalDir;

    std::vector<std::string> fCalibrationFile;
    std::vector<std::string> fRawCalibrationFile;

};












#endif
