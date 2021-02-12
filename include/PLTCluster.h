#ifndef GUARD_PLTCluster_h
#define GUARD_PLTCluster_h

#include "PLTHit.h"

#include <iostream>
#include <vector>

class PLTCluster
{
  public:
    PLTCluster ();
    ~PLTCluster ();

    void AddHit (PLTHit*);
    float Charge ();
    size_t NHits ();
    PLTHit* Hit (size_t const);
    PLTHit* SeedHit ();


    int PX ();
    int PY ();
    int PZ ();
    std::pair<int, int> PCenter ();
    int Channel ();
    int ROC ();

    // Local w.r.t. center of diamond
    float LX ();
    float LY ();
    std::pair<float, float> LCenter ();

    // Telescope cords
    float TX ();
    float TY ();
    float TZ ();
    std::pair<float, float> TCenter ();

    // Global cords
    float GX ();
    float GY ();
    float GZ ();
    std::pair<float, float> GCenter ();

    // Cluster center
    std::pair<float, float> CenterOfMass(std::string);
    std::pair<float, float> LCenterOfMass() { return CenterOfMass("local"); }
    std::pair<float, float> TCenterOfMass() { return CenterOfMass("telescope"); }
    std::pair<float, float> GCenterOfMass() { return CenterOfMass("global"); }


  private:
    std::vector<PLTHit*> fHits;  // The seed hit needs to be 0 in this vector

};




#endif
