#ifndef GUARD_PLTTrack_h
#define GUARD_PLTTrack_h

#include <vector>
#include <iostream>
#include <math.h>

#include "TGraph.h"
#include "TGraphErrors.h"

#include "TF1.h"


#include "PLTCluster.h"
#include "PLTAlignment.h"
#include "PLTPlane.h"
#include "PLTU.h"


class PLTTrack
{
  public:
    PLTTrack ();
    ~PLTTrack ();

    void AddCluster (PLTCluster*);
    int  MakeTrack (PLTAlignment&, int);

    size_t NClusters() { return fClusters.size(); }
    size_t NHits();

    PLTCluster * Cluster (size_t const i ) { return fClusters[i]; }
    float LResidualX (size_t const i) { return fLResidualX[i]; }
    float LResidualY (size_t const i) { return fLResidualY[i]; }
    std::pair<float, float> GetResiduals (PLTCluster&, PLTAlignment &);
    std::pair<float, float> GetResiduals(size_t const i) { return std::make_pair(fLResidualX[i], fLResidualY[i]); }

    bool IsFiducial (PLTPlane*, PLTAlignment&, PLTPlane::FiducialRegion);
    bool IsFiducial (int, int, PLTAlignment&, PLTPlane::FiducialRegion);

    float TX (float);
    float TY (float);

    std::pair<float, float> GXYatGZ (float, PLTAlignment&);

    float D2 ();

    float Chi2(){return fChi2;}

    float Chi2X(){return fChi2X;}

    float Chi2Y(){return fChi2Y;}

    float ExtrapolateX(float Z){return (Z * tan(fAngleRadX)) + fOffsetX;}
    float ExtrapolateY(float Z){return (Z * tan(fAngleRadY)) + fOffsetY;}

    float InterPolateX(float z) { return z * fSlopeX + fOffsetX; }
    float InterPolateY(float z) { return z * fSlopeY + fOffsetY; }

  private:
    std::vector<PLTCluster*> fClusters;

  public:

    // Vector in *telescope* and *global* coords
    float fTVX;
    float fTVY;
    float fTVZ;

    float fGVX;
    float fGVY;
    float fGVZ;

    // Origin in *telescope* and *global* coords as defined by ROC-0
    float fTOX;
    float fTOY;
    float fTOZ;

    float fGOX;
    float fGOY;
    float fGOZ;

    // Track Fit result for slope (use to get telescope coordinates for any point in z-direction)
    float fAngleX, fAngleY, fOffsetX, fOffsetY, fAngleRadX, fAngleRadY;
    float fSlopeX, fSlopeY;

    // Where the track passes through the X=0(=0), Y=0(=1), and Z=0 planes
    // Three corrds just because that's easy enough
    float fPlaner[3][3];

    // Residuals for each roc in X and Y in terms of pixels
    std::vector<float> fLResidualX;
    std::vector<float> fLResidualY;

    float fD2, fChi2, fChi2X, fChi2Y;

    static bool const DEBUG = false;

};

#endif
