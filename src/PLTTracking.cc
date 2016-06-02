#include "PLTTracking.h"

// Needed these for the track-creation loop

// Types to hold vector-of-clusters (Vc) and vector-of-vector-of-clusters (Vvc)
typedef std::vector<PLTCluster*> Vc;
typedef std::vector<Vc> Vvc;

// Use a vector of iterators
// which iterate over the individual vector<clusters>s.
struct Digits {
    Vc::const_iterator begin;
    Vc::const_iterator end;
    Vc::const_iterator me;
};
typedef std::vector<Digits> Vd;


PLTTracking::PLTTracking (int nplanes) : fNPlanes(nplanes)
{
  SetAllPlanes();
}

PLTTracking::~PLTTracking ()
{
}


void PLTTracking::SetTrackingAlignment (PLTAlignment* Alignment)
{
  if (!Alignment) {
    std::cerr << "ERROR: PLTTracking::SetTrackingAlignment(): Alignment object does not exist: is 0x0" << std::endl;
    exit(1);
  }

  fAlignment = Alignment;
  return;
}


void PLTTracking::SetTrackingAlgorithm (TrackingAlgorithm const Algorithm)
{
  fTrackingAlgorithm = Algorithm;
  return;
}


int PLTTracking::GetTrackingAlgorithm ()
{
  return fTrackingAlgorithm;
}


void PLTTracking::SetAllPlanes(){

  fUsePlanesForTracking.resize(0);
  for (int i=0;i!=fNPlanes;i++) {
    fUsePlanesForTracking.push_back(2);
  }

  fDoSinglePlaneEfficiency = false;
}


void PLTTracking::SetPlaneUnderTest(int put){

    fDoSinglePlaneEfficiency = true;

    // The default is 222222 -> require at least one hit in all planes
    // (there and an additional condition before calling tracking that
    //  requires exactly 6 hits)

    // For single-plane efficiency we want (for example): 333033
    for (uint8_t i=0;i!=fUsePlanesForTracking.size();i++){
      if (i==put){
        fUsePlanesForTracking[i] = 0;
      }
      else {
        fUsePlanesForTracking[i] = 3;
      }
    }
}


void PLTTracking::SetPlaneUnderTestSandwich( int put){

    fDoSinglePlaneEfficiency = true;

    // The default is 222222 -> require at least one hit in all planes
    // (there and an additional condition before calling tracking that
    //  requires exactly 6 hits)

    // For sandwich measurement (to get the residual) we want (for example):
    // 303000
    if (put==0){
      fUsePlanesForTracking[0] = 0;
      fUsePlanesForTracking[1] = 3;
      fUsePlanesForTracking[2] = 3;
      fUsePlanesForTracking[3] = 0;
      fUsePlanesForTracking[4] = 0;
      fUsePlanesForTracking[5] = 0;
    }
    else if (put==1){
      fUsePlanesForTracking[0] = 3;
      fUsePlanesForTracking[1] = 0;
      fUsePlanesForTracking[2] = 3;
      fUsePlanesForTracking[3] = 0;
      fUsePlanesForTracking[4] = 0;
      fUsePlanesForTracking[5] = 0;
    }
    else if (put==2){
      fUsePlanesForTracking[0] = 0;
      fUsePlanesForTracking[1] = 3;
      fUsePlanesForTracking[2] = 0;
      fUsePlanesForTracking[3] = 3;
      fUsePlanesForTracking[4] = 0;
      fUsePlanesForTracking[5] = 0;
    }
    else if (put==3){
      fUsePlanesForTracking[0] = 0;
      fUsePlanesForTracking[1] = 0;
      fUsePlanesForTracking[2] = 3;
      fUsePlanesForTracking[3] = 0;
      fUsePlanesForTracking[4] = 3;
      fUsePlanesForTracking[5] = 0;
    }
    else if (put==4){
      fUsePlanesForTracking[0] = 0;
      fUsePlanesForTracking[1] = 0;
      fUsePlanesForTracking[2] = 0;
      fUsePlanesForTracking[3] = 3;
      fUsePlanesForTracking[4] = 0;
      fUsePlanesForTracking[5] = 3;
    }
    else if (put==5){
      fUsePlanesForTracking[0] = 0;
      fUsePlanesForTracking[1] = 0;
      fUsePlanesForTracking[2] = 0;
      fUsePlanesForTracking[3] = 3;
      fUsePlanesForTracking[4] = 3;
      fUsePlanesForTracking[5] = 0;
    }



}


void PLTTracking::RunTracking (PLTTelescope& Telescope)
{
  switch (fTrackingAlgorithm) {
    case kTrackingAlgorithm_NoTracking:
      break;
    case kTrackingAlgorithm_01to2_All:
      TrackFinder_01to2_All(Telescope);
      break;
    case kTrackingAlgorithm_6PlanesHit:
      TrackFinder_AllPlanesHit(Telescope);
      break;
    case kTrackingAlgorithm_ETH:
      TrackFinder_ETH(Telescope);
      break;
    default:
      std::cerr << "ERROR: PLTTracking::RunTracking() has no idea what tracking algorithm you want to use" << std::endl;
      throw;
  }

  return;
}


bool PLTTracking::CompareTrackD2 (PLTTrack* lhs, PLTTrack* rhs)
{
  return lhs->D2() < rhs->D2();
}



void PLTTracking::TrackFinder_01to2_All (PLTTelescope& Telescope)
{
  // Find tracks in this telescope which are more or less parallel to beam
  // Beamspot assumed to be at 0,0,0.

  // Check that tracks haven't already been filled
  if (Telescope.NTracks() != 0) {
    std::cerr << "ERROR: It looks like tracks have already been filled here: PLTTelescope::TrackFinderParallelTracks()" << std::endl;
    return;
  }

  // If there isn't a hit in each plane skip it for now
  if (Telescope.HitPlaneBits() != 0x7) {
    return;
  }

  std::vector<int> UsedHits;

  // Shorthand for each plane
  PLTPlane* P0 = Telescope.Plane(0);
  PLTPlane* P1 = Telescope.Plane(1);
  PLTPlane* P2 = Telescope.Plane(2);

  if (P0->NClusters() == 0 || P1->NClusters() == 0 || P2->NClusters() == 0) {
    return;
  }

  // Vector to keep track of tracks that we're interested in
  std::vector<PLTTrack*> MyTracks;

  // Start seeding with clusters in the 0th plane
  for (size_t iCL0 = 0; iCL0 != P0->NClusters(); ++iCL0) {


    // Loop over two other planes/clusters.  Pick a cluster in 1st plane,
    // and see if any patch behind in the 2nd plane.  One can rank these by Chi2 and pick best
    for (size_t iCL1 = 0; iCL1 != P1->NClusters(); ++iCL1) {

      PLTTrack Track01;
      Track01.AddCluster(P0->Cluster(iCL0));
      Track01.AddCluster(P1->Cluster(iCL1));
      Track01.MakeTrack(*fAlignment, Telescope.NPlanes() );

      float const ProjectionX2 = Track01.TX(P2->TZ());
      float const ProjectionY2 = Track01.TY(P2->TZ());

      for (size_t iCL2 = 0; iCL2 != P2->NClusters(); ++iCL2) {
        float const X = P2->Cluster(iCL2)->TX();
        float const Y = P2->Cluster(iCL2)->TY();

        float const ResidualX = ProjectionX2 - X;
        float const ResidualY = ProjectionY2 - Y;

        float const Distance = sqrt(ResidualX*ResidualX + ResidualY*ResidualY);

        if (DEBUG) {
          printf("Distance: %12.3E\n", Distance);
        }


        // If it's not too far off, keep it!
        if (Distance < 0.2000) {
          // Keep as possible track..
          PLTTrack* Track012 = new PLTTrack();
          Track012->AddCluster(P0->Cluster(iCL0));
          Track012->AddCluster(P1->Cluster(iCL1));
          Track012->AddCluster(P2->Cluster(iCL2));
          Track012->MakeTrack(*fAlignment, Telescope.NPlanes());
          MyTracks.push_back(Track012);
        }

      }
    }

  }

  int const NTracksBefore = (int) MyTracks.size();

  // Grab the best tracks first and don't let there be overlap..
  SortOutTracksNoOverlapBestD2(MyTracks);

  if (DEBUG) {
    printf("Found NTracks possible: %4i   Kept NTracks: %4i\n", NTracksBefore, (int) MyTracks.size());
  }

  for (size_t i = 0; i != MyTracks.size(); ++i) {
    Telescope.AddTrack(MyTracks[i]);
  }

  return;
}


void PLTTracking::TrackFinder_AllPlanesHit (PLTTelescope& Telescope)
{
    /** Find tracks in this telescope which are more or less parallel to beam
        Beamspot assumed to be at 0,0,0.

        Check that tracks haven't already been filled */
    if (Telescope.NTracks() != 0){
        std::cerr << "ERROR: It looks like tracks have already been filled here: PLTTelescope::TrackFinderParallelTracks()" << std::endl;
        return;
    }

    /** Check if all the planes with mandatory clusters also have a cluster */
    for (uint8_t iPlane=0; iPlane < Telescope.NPlanes(); iPlane++){
        if ( (fUsePlanesForTracking[iPlane]==2) && (Telescope.Plane(iPlane)->NClusters()==0) ){
            return;
        }
    }

    // Check if all the planes which require exactly one cluster have that
    for (uint8_t iPlane=0; iPlane < Telescope.NPlanes(); iPlane++){
        if ((fUsePlanesForTracking[iPlane]==3) &&
                (Telescope.Plane(iPlane)->NClusters() != 1 ))
        {
            return;
        }
    }

  // Put all the clusters we actually use for
  // tracking into a matrix.
  // Need to have fUsePlanesForTracking of either 1 or 2 or 3
  //  and > 0 hits
  std::vector< std::vector< PLTCluster* > > ClustersForTracking;
  for (uint8_t iPlane=0; iPlane < Telescope.NPlanes(); iPlane++){
    if ( (fUsePlanesForTracking[iPlane]>0)
      && (Telescope.Plane(iPlane)->NClusters()>0)){cout << "/n bla si /n" << endl;

        std::vector< PLTCluster* > VClusters;
        for (size_t iCluster = 0;
            iCluster != Telescope.Plane(iPlane)->NClusters();
            ++iCluster) {
              VClusters.push_back( Telescope.Plane(iPlane)->Cluster(iCluster) );
        } // end: loop over clusters

        ClustersForTracking.push_back( VClusters );

    } // end: use plane and plane has clusters
  } // end: loop over planes

  // Vector to keep track of tracks that we're interested in
  std::vector<PLTTrack*> MyTracks;

  // Code adapted from:
  // http://stackoverflow.com/questions/5279051/how-can-i-create-cartesian-product-of-vector-of-vectors
  // Idea is to have a vector of "Digit" objects that store the information to
  // loop over a vector<cluster*>.
  // Start all of the iterators at the beginning.
  Vd vd;
  for(Vvc::const_iterator it = ClustersForTracking.begin();
    it != ClustersForTracking.end();
    ++it) {
    Digits d = {(*it).begin(), (*it).end(), (*it).begin()};
    vd.push_back(d);
  } // end of initializing the digits

  // Make sure we have at least two planes with hits
  // otherwise it will be hard to make tracks
  if (vd.size() < 2)
    return;

  // Actual track creation
  bool keepRunning = true;
  while(keepRunning) {

    PLTTrack* Track = new PLTTrack();
    // Construct the firsr track track by accessing
    // the iterators in the Vd object
    for(Vd::const_iterator it = vd.begin(); it != vd.end(); it++)
        Track->AddCluster(*(it->me));

    Track->MakeTrack(*fAlignment, Telescope.NPlanes() );
    MyTracks.push_back(Track);


    // Increment the Digits
    for(Vd::iterator it = vd.begin(); ; ) {
        // Start with the left one
        ++(it->me);
        // If it hits its end:
        if(it->me == it->end) {
          // test if we are at the right end
          if(it+1 == vd.end()) {
            keepRunning = false;
            break;
          // Otherwise cascade (set to zero, increment next one on the right)
          } else {
              it->me = it->begin;
              ++it;
            }
        // Normal increase - no cascading necessary
        } else {
            break;
          }
    } // end of incrementing the digits
  } // end track creation

  int const NTracksBefore = (int) MyTracks.size();

  // Grab the best tracks first and don't let there be overlap..
  SortOutTracksNoOverlapBestD2(MyTracks);

  if (DEBUG) {
    printf("Found NTracks possible: %4i   Kept NTracks: %4i\n", NTracksBefore, (int) MyTracks.size());
  }

  for (size_t i = 0; i != MyTracks.size(); ++i) {
    Telescope.AddTrack(MyTracks[i]);
  }

  return;
}

void PLTTracking::TrackFinder_ETH (PLTTelescope& Telescope)
{
  // Find tracks in this telescope which are more or less parallel to beam
  // Beamspot assumed to be at 0,0,0.

  // Check that tracks haven't already been filled
  if (Telescope.NTracks() != 0) {
    std::cerr << "ERROR: It looks like tracks have already been filled here: PLTTelescope::TrackFinderParallelTracks()" << std::endl;
    return;
  }

  // Check if the first four planes with mandatory clusters also have a cluster
  for (uint8_t iPlane=0; iPlane < 4; iPlane++){
    if ( (fUsePlanesForTracking[iPlane]==2)
      && (Telescope.Plane(iPlane)->NClusters()==0)){
        return;
    }
  }
//  /** just for hte moment only one hit */
//    for (uint8_t iPlane=0; iPlane < 4; iPlane++){
//    if ( (fUsePlanesForTracking[iPlane]==2)
//      && (Telescope.Plane(iPlane)->NHits() !=1)){
//        return;
//    }
//  }

///** take out weird tracks */
//    if ( (fUsePlanesForTracking[0]==2) && (Telescope.Plane(0)->Hit(0)->Column() - Telescope.Plane(1)->Hit(0)->Column() > 20) ){
//        return;
//    }

  // Check if all the planes which require exactly one cluster have that
  for (uint8_t iPlane=0; iPlane < 4; iPlane++){
    if ((fUsePlanesForTracking[iPlane]==3) &&
        (Telescope.Plane(iPlane)->NClusters() != 1 )){
        return;
    }
  }

  // Put all the clusters we actually use for
  // tracking into a matrix.
  // Need to have fUsePlanesForTracking of either 1 or 2
  //  and > 0 hits
  std::vector< std::vector< PLTCluster* > > ClustersForTracking;
  for (uint8_t iPlane=0; iPlane < 4; iPlane++){
    if ( (fUsePlanesForTracking[iPlane]>0)
      && (Telescope.Plane(iPlane)->NClusters()>0)){

        std::vector< PLTCluster* > VClusters;
        for (size_t iCluster = 0;
            iCluster != Telescope.Plane(iPlane)->NClusters();
            ++iCluster) {
              VClusters.push_back( Telescope.Plane(iPlane)->Cluster(iCluster) );
        } // end: loop over clusters

        ClustersForTracking.push_back( VClusters );

    } // end: use plane and plane has clusters
  } // end: loop over planes

  // Vector to keep track of tracks that we're interested in
  std::vector<PLTTrack*> MyTracks;

  // Code adapted from:
  // http://stackoverflow.com/questions/5279051/how-can-i-create-cartesian-product-of-vector-of-vectors
  // Idea is to have a vector of "Digit" objects that store the information to
  // loop over a vector<cluster*>.
  // Start all of the iterators at the beginning.
  Vd vd;
  for(Vvc::const_iterator it = ClustersForTracking.begin();
    it != ClustersForTracking.end();
    ++it) {
    Digits d = {(*it).begin(), (*it).end(), (*it).begin()};
    vd.push_back(d);
  } // end of initializing the digits

  // Make sure we have at least two planes with hits
  // otherwise it will be hard to make tracks
  if (vd.size() < 2)
    return;

  // Actual track creation
  bool keepRunning = true;
  while(keepRunning) {

    PLTTrack* Track = new PLTTrack();
    // Construct the firsr track track by accessing
    // the iterators in the Vd object
    for(Vd::const_iterator it = vd.begin(); it != vd.end(); it++)
        Track->AddCluster(*(it->me));

    Track->MakeTrack(*fAlignment, 4 );
    MyTracks.push_back(Track);


    // Increment the Digits
    for(Vd::iterator it = vd.begin(); ; ) {
        // Start with the left one
        ++(it->me);
        // If it hits its end:
        if(it->me == it->end) {
          // test if we are at the right end
          if(it+1 == vd.end()) {
            keepRunning = false;
            break;
          // Otherwise cascade (set to zero, increment next one on the right)
          } else {
              it->me = it->begin;
              ++it;
            }
        // Normal increase - no cascading necessary
        } else {
            break;
          }
    } // end of incrementing the digits
  } // end track creation

  int const NTracksBefore = (int) MyTracks.size();

  // Grab the best tracks first and don't let there be overlap..
  SortOutTracksNoOverlapBestD2(MyTracks);

  if (DEBUG) {
    printf("Found NTracks possible: %4i   Kept NTracks: %4i\n", NTracksBefore, (int) MyTracks.size());
  }

  for (size_t i = 0; i != MyTracks.size(); ++i) {
    Telescope.AddTrack(MyTracks[i]);
  }
//  std::cout << "TRACK:" << std::endl;
//  for (uint8_t i = 0; i != Telescope.Track(0)->NClusters(); i++ ){
//    if (Telescope.Track(0)->Cluster(0)->TX() - Telescope.Track(0)->Cluster(1)->TX() > 0.2){
//        std::cout.precision(3);
//        std::cout << Telescope.Track(0)->Cluster(i)->LX() << "\t" << Telescope.Track(0)->Cluster(i)->TX() << "\t" << Telescope.Track(0)->Cluster(i)->Hit(0)->Column();
//        std::cout << "\t" << Telescope.Track(0)->Cluster(i)->Hit(0)->ROC() << std::endl;
//        std::cout << Telescope.Track(0)->Cluster(i)->LY() << "\t" << Telescope.Track(0)->Cluster(i)->TY() << "\t" << Telescope.Track(0)->Cluster(i)->Hit(0)->Row() << std::endl;
//    }
//  }
//  std::cout << std::endl;
  return;
}







void PLTTracking::SortOutTracksNoOverlapBestD2 (std::vector<PLTTrack*>& MyTracks)
{
  // Idea of this function is to start with the tracks with the best test-stat
  // and grab those tracks first..  then for the remaining tracks only keep them
  // if they have unique clusters.

  // Order tracks with lowest test statistic first
  std::sort(MyTracks.begin(), MyTracks.end(), &PLTTracking::CompareTrackD2);

  // Containers for tracks and clsters
  std::vector<PLTTrack*> UsedTracks;
  std::vector<PLTTrack*> SkippedTracks;
  std::set<PLTCluster*> UsedClusters;

  // Loop over all tracks so far
  for (std::vector<PLTTrack*>::iterator Track = MyTracks.begin(); Track != MyTracks.end(); ++Track) {
    int UsedClusterFlag = false;
    for (size_t icluster = 0; icluster != (*Track)->NClusters(); ++icluster) {
      if (UsedClusters.count((*Track)->Cluster(icluster))) {
        UsedClusterFlag = true;
        break;
      }
    }
    if (!UsedClusterFlag) {

      for (size_t icluster = 0; icluster != (*Track)->NClusters(); ++icluster) {

        // If it's a track with all unused clusters so far let's keep it
        UsedClusters.insert((*Track)->Cluster(icluster));
      }
      UsedTracks.push_back(*Track);
    } else {
      // This track has a cluster which has already been used, remember it
      // so we can delete it later
      SkippedTracks.push_back(*Track);
    }
  }

  // Replace input vector of tracks with accepted tracks
  MyTracks = UsedTracks;

  // Delete the unused tracks
  for (size_t i = 0; i != SkippedTracks.size(); ++i) {
    delete SkippedTracks[i];
  }

  return;
}

