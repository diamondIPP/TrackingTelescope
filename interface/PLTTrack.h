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

    size_t NClusters ();
    size_t NHits ();

    PLTCluster* Cluster (size_t const);

    float LResidualX (size_t const);
    float LResidualY (size_t const);
    std::pair<float, float> LResiduals (PLTCluster&, PLTAlignment&);

    bool IsFiducial (PLTPlane*, PLTAlignment&, PLTPlane::FiducialRegion);
    bool IsFiducial (int const, int const, PLTAlignment&, PLTPlane::FiducialRegion);

    float TX (float const);
    float TY (float const);

    std::pair<float, float> GXYatGZ (float const, PLTAlignment&);

    float D2 ();

    float Chi2(){return fChi2;}

    float Chi2X(){return fChi2X;}

    float Chi2Y(){return fChi2Y;}

    float ExtrapolateX(float Z){return (Z * tan(fAngleRadX)) + fOffsetX;}
    float ExtrapolateY(float Z){return (Z * tan(fAngleRadY)) + fOffsetY;}

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
