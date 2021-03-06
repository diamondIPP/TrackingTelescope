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
    PSIFileReader(bool track_only_telescope);
    virtual ~PSIFileReader() = default;

    virtual bool OpenFile () = 0;
    virtual void ResetFile () = 0 ;

    void Clear ();

    void ReadPixelMask (std::string const);
    void AddToPixelMask( int, int, int, int);
    bool IsPixelMasked (int const);

    void DrawTracksAndHits (std::string const);

    virtual int GetNextEvent () = 0;
    virtual unsigned GetEntries() = 0;
    virtual void CloseFile() = 0;

    size_t NHits ();
    PLTHit* Hit (size_t);

    PLTGainCal * GetGainCal () { return &fGainCal; }
    PLTAlignment * GetAlignment() { return &fAlignment; }
    const std::set<int> * GetPixelMask(){ return &fPixelMask; }

protected:

    std::set<int> fPixelMask;
    std::vector<PLTHit*> fHits;

    PSIGainInterpolator fGainInterpolator;
    PLTGainCal fGainCal;
    PLTAlignment fAlignment;

    std::string fBinaryFileName;

    std::map<int, PLTPlane> fPlaneMap;

    std::string fBaseCalDir;

    std::vector<std::string> fCalibrationFile;
    std::vector<std::string> fRawCalibrationFile;

};












#endif
