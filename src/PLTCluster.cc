#include "PLTCluster.h"


PLTCluster::PLTCluster ()
{
  // Make it real good
}


PLTCluster::~PLTCluster ()
{
  // Makeun me
}



void PLTCluster::AddHit (PLTHit* Hit)
{
  // Add a hit
  fHits.push_back(Hit);
  return;
}


float PLTCluster::Charge ()
{
  // Compute charge of this cluster
  float Sum = 0;
  for (std::vector<PLTHit*>::iterator it = fHits.begin(); it != fHits.end(); ++it) {
    Sum += (*it)->Charge();
  }

  return Sum;
}


size_t PLTCluster::NHits ()
{
  return fHits.size();
}


PLTHit* PLTCluster::Hit (size_t const i)
{
  return fHits[i];
}

PLTHit* PLTCluster::SeedHit ()
{
  return fHits[0];
}


int PLTCluster::PX ()
{
  return PCenter().first;
}


int PLTCluster::PY ()
{
  return PCenter().second;
}


int PLTCluster::PZ ()
{
  return SeedHit()->ROC();
}


std::pair<int, int> PLTCluster::PCenter ()
{
  return std::make_pair<int, int>(SeedHit()->Column(), SeedHit()->Row());
}



int PLTCluster::Channel ()
{
  return SeedHit()->Channel();
}


int PLTCluster::ROC ()
{
  return SeedHit()->ROC();
}




float PLTCluster::LX ()
{
  return LCenter().first;
}


float PLTCluster::LY ()
{
  return LCenter().second;
}



std::pair<float, float> PLTCluster::LCenter ()
{
  return LCenterOfMass();
  //return std::make_pair<float, float>(SeedHit()->LX(), SeedHit()->LY());
}


float PLTCluster::TX ()
{
  return TCenter().first;
}


float PLTCluster::TY ()
{
  return TCenter().second;
}


float PLTCluster::TZ ()
{
  return SeedHit()->TZ();
}


std::pair<float, float> PLTCluster::TCenter ()
{

  return TCenterOfMass();
}

float PLTCluster::GX ()
{
  return GCenter().first;
}


float PLTCluster::GY ()
{
  return GCenter().second;
}


float PLTCluster::GZ ()
{
  return SeedHit()->GZ();
}


std::pair<float, float> PLTCluster::GCenter ()
{
  return GCenterOfMass();
  //return std::make_pair<float, float>(SeedHit()->GX(), SeedHit()->GY());
}

std::pair<float, float> PLTCluster::CenterOfMass(std::string type) {

  /** Return the coordinates based on a charge weighted average of pixel hits

     X, Y, and Total Charge */
  float X(0.0), Y(0.0);
  float ChargeSum(0.0);
  bool FoundZeroCharge = false;

  /** We should never have negative or zero charges and cannot handle them -> just take the unweighted average in that case */
  for (auto fHit: fHits)
    if (fHit->Charge() == -9999 or fHit->Charge() < 0)
      FoundZeroCharge = true;

  // Loop over each hit in the cluster
  for (auto &fHit : fHits) {
    float iCharge = (FoundZeroCharge ? 1 : fHit->Charge());
    float iX, iY;
    if (type.find("local") != std::string::npos){
      iX = fHit->LX();
      iY = fHit->LY();
    }
    else if ((type.find("global") != std::string::npos)){
      iX = fHit->GX();
      iY = fHit->GY();
    }
    else {
      iX = fHit->TX();
      iY = fHit->TY();
    }

    X += iX * iCharge;
    Y += iY * iCharge;
    ChargeSum += fHit->Charge();
  }

  /** If charge sum is zero or less or we have a zero charge return average */
  if (ChargeSum <= 0.0 or FoundZeroCharge)
    return std::make_pair<float, float>(X / (float) NHits(), Y / (float) NHits());

  return std::make_pair<float, float>(X / ChargeSum, Y / ChargeSum);
}
