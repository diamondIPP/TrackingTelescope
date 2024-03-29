#include "PLTAnalysis.h"
#include "Utils.h"

using namespace std;

PLTAnalysis::PLTAnalysis(string const & inFileName, TFile * Out_f,  TString const & runNumber, uint8_t const TelescopeID, bool TrackOnlyTelescope, uint64_t max_event_nr):
    Action(inFileName, runNumber),
    telescopeID(TelescopeID),
    now1(clock()), now2(clock()), loop(0), startProg(0), endProg(0), allProg(0), averTime(0),
    TimeWidth(20000), StartTime(0), NGraphPoints(0),
    PHThreshold(3e5), is_root_file_(IsROOTFile(inFileName)), trackOnlyTelescope(TrackOnlyTelescope)
{
    out_f = Out_f;
    /** set up root */
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;
    /** single plane studies */
    SinglePlaneStudies();
    /** init file reader */
    FR = InitFileReader();
    if (is_root_file_) nEntries = ((PSIRootFileReader*) FR)->fTree->GetEntries();
    stopAt = max_event_nr ? max_event_nr : nEntries;
    /** apply masking */
    FR->ReadPixelMask(GetMaskingFilename());
    /** init histos */
    Histos = new RootItems(run_number_);
    DiaZ = getDiaZPositions();
    cout << "Output directory: " << Histos->getOutDir() << endl;
    /** init file writer */
    if (UseFileWriter())
      FW = new FileWriterTracking(in_file_name_, FR);
    PBar = new tel::ProgressBar(stopAt - 1);
}

PLTAnalysis::~PLTAnalysis()
{
  if (UseFileWriter()) {
    FW->saveTree();
    delete FR;
  }
    delete FW;
//    delete Histos; // This causes it to crash for some unknown reason...
}


/** ============================
 EVENT LOOP
 =================================*/
 void PLTAnalysis::EventLoop(){

    getTime(now1, startProg);
    now1 = clock();
//    cout << "stopAt = " << stopAt << endl;
//        stopAt = 1e5;
    for (uint32_t ievent = 0; FR->GetNextEvent() >= 0; ++ievent) {
        if (ievent > stopAt) break;
        ThisTime = ievent;

        PBar->update(ievent);
        /** file writer */
        if (is_root_file_) { WriteTrackingTree(); }

        /** fill coincidence map */
        Histos->CoincidenceMap()->Fill(FR->HitPlaneBits() );

        /** make average pulseheight maps*/
        if (ThisTime - (StartTime + NGraphPoints * TimeWidth) > TimeWidth)
            MakeAvgPH();

        /** draw tracks if there is more than one hit*/
        DrawTracks();

        /** loop over the planes */
        for (uint8_t iplane = 0; iplane != FR->NPlanes(); ++iplane) {

            /** quick alignment check */
//            bool aligned = FR->NTracks() > 0 ? FR->NTracks() and FR->Plane(iplane)->NHits() : true;  //set automatically true if there is no track (don't care about these events)
//            FW->setAligned(iplane, FW->lastIsAligned(iplane) or aligned);  //only set misaligned mark if this and the last event are not aligned
//            FW->setOldAligned(iplane, aligned);

            PLTPlane * Plane = FR->Plane(iplane);
            /** Check that the each hit belongs to only one cluster type*/ //todo: DA: comentar
//    	    Plane->CheckDoubleClassification();
            /** fill cluster histo */
            Histos->nClusters()[Plane->ROC()]->Fill(Plane->NClusters());

            for (size_t icluster = 0; icluster != Plane->NClusters(); ++icluster) {
                PLTCluster* Cluster = Plane->Cluster(icluster);

                /** ignore charges above a certain threshold */
                if (Cluster->Charge() > PHThreshold) continue;

                /** fill pulse heigh histos */
                FillPHHistos(iplane, Cluster);

                /** fill hits per cluster histo */
                Histos->nHitsPerCluster()[Cluster->ROC()]->Fill(Cluster->NHits());

                /** fill high and low occupancy */
                FillOccupancyHiLo(Cluster);
            }
            /** fill occupancy histo */
            FillOccupancy(Plane);
        }


//        cout << "Number of Tracks: " << FR->NTracks() << endl;
        if (UseFileWriter() && FR->NTracks() == 1 && ((FR->Track(0)->NClusters() == Histos->NRoc()  && !trackOnlyTelescope) || (FR->Track(0)->NClusters() >= 4  && trackOnlyTelescope))){
//		if ((telescopeID == 7 || telescopeID == 8 || telescopeID == 9 || telescopeID >= 10) && FR->NTracks() == 1 && FR->Track(0)->NClusters() == Histos->NRoc() ){

		  do_slope = true;
		  for (uint8_t i_rocs(0); i_rocs != FR->Track(0)->NClusters(); i_rocs++)
		    if (FR->Track(0)->Cluster(i_rocs)->Charge() > PHThreshold){
		      do_slope = false;
		      break;
		    }

		  if (do_slope) {

		    PLTTrack * Track = FR->Track(0);
		    //				for (uint8_t i=0; i != FR->Signal().size(); i++)
		    //                    std::cout << FR->SignalDiamond(i) << " ";
		    //                std::cout << std::endl;

		    //                if (Track->NHits() == 4){
		    //                    for (uint16_t i = 0; i < Track->NClusters(); i++){
		    //                        cout << "Plane: " << i << ":\t" << setprecision(2) << setw(4) << Track->Cluster(i)->TX()*100 << "\t" << Track->Cluster(i)->TY()*100;
		    //                        cout << "\t" << Track->Cluster(i)->Hit(0)->Column() << "\t"<< Track->Cluster(i)->Hit(0)->Row() << endl;
		    //                    }
		    //                    cout << Track->fChi2X << " " << Track->fChi2Y << " " << Track->fChi2 << endl;
		    ////                    float y1 = Track->Cluster(2)->TY() - Track->Cluster(1)->TY();
		    ////                    float y2 = Track->Cluster(1)->TY() - Track->Cluster(0)->TY();
		    ////                    cout << y1*100 << " " << y2*100 << endl;
		    ////                    cout << (Track->Cluster(2)->TX() - Track->Cluster(1)->TX()) / (Track->Cluster(1)->TX() - Track->Cluster(0)->TX()) * 2.032 << " ";
		    ////                    cout << (Track->Cluster(2)->TY() - Track->Cluster(1)->TY()) / (Track->Cluster(1)->TY() - Track->Cluster(0)->TY()) * 2.032 << endl<<endl;
		    //
		    //                }

		    /** fill chi2 histos */
		    Histos->Chi2()->Fill(Track->Chi2());
		    Histos->Chi2X()->Fill(Track->Chi2X());
		    Histos->Chi2Y()->Fill(Track->Chi2Y());

		    /** fill slope histos */
		    Histos->TrackSlopeX()->Fill(Track->fAngleX);
		    Histos->TrackSlopeY()->Fill(Track->fAngleY);

		    //                if (ievent < 100){
		    //                    for (uint8_t iSig = 0; iSig != Track->NClusters(); iSig++)
		    //                        std::cout<< Track->Cluster(iSig)->TX() << " " << Track->Cluster(iSig)->TY() << " " << Track->Cluster(iSig)->TZ() << std::endl;
		    //                        std::cout << Track->fChi2X << " " << Track->fChi2Y << " " << Track->fChi2 << std::endl;
		    //                        std::cout << Track->fAngleRadX << " " << Track->fOffsetX << Track->fAngleRadY << " " << Track->fOffsetY << std::endl;
		    //                    std::cout << std::endl;
		    //                }

		    /** fill signal histos */
		    if (FillSignalHistos()) {
//		    if (telescopeID == 9 || telescopeID == 8 || telescopeID ==7) {
		      if (ievent > 0 && FW->InTree()->GetBranch(GetSignalBranchName())){
                        for (uint8_t iSig = 0; iSig != Histos->NSig(); iSig++){
              float dia1z = tel::Config::dia_z_pos_.at(0);
              float dia2z = tel::Config::dia_z_pos_.at(1);
			  if (iSig < 2)
			    Histos->SignalDisto()[iSig]->Fill(Track->ExtrapolateX(dia1z), Track->ExtrapolateY(dia1z), FR->SignalDiamond(iSig) );
			  else
			    Histos->SignalDisto()[iSig]->Fill(Track->ExtrapolateX(dia2z), Track->ExtrapolateY(dia2z), FR->SignalDiamond(iSig) );
                        }
		      }
		    }

		    /** loop over the clusters */
		    for (size_t icluster = 0; icluster < Track->NClusters(); icluster++){

					/** get the ROC in of the cluster and fill the corresponding residual */
		      uint8_t ROC = Track->Cluster(icluster)->ROC();
		      PLTCluster * Cluster = Track->Cluster(icluster);

		      /** fill residuals */
		      Histos->Residual()[ROC]->Fill(Track->LResidualX(ROC), Track->LResidualY(ROC)); // dX vs dY
					Histos->ResidualXdY()[ROC]->Fill(Cluster->LX(), Track->LResidualY(ROC));// X vs dY
					Histos->ResidualYdX()[ROC]->Fill(Cluster->LY(), Track->LResidualX(ROC)); // Y vs dX

					/** ignore events above a certain threshold */
					if (Cluster->Charge() > PHThreshold) continue;
					//printf("High Charge: %13.3E\n", Cluster->Charge());

					/** fill the offline pulse heights (Track6+|Slope| < 0.01 in x and y ) */
					FillOfflinePH(Track, Cluster);
		    }
		  }
		}


    } /** END OF EVENT LOOP */
    /** add the last point to the average pulse height graph */
    MakeAvgPH();

    cout << endl;
    getTime(now1, loop);
    now1 = clock();
 }
/** ============================
 AFTER LOOP -> FINISH
 =================================*/
 void PLTAnalysis::FinishAnalysis(){

    out_f->cd();

    Histos->SaveAllHistos();

    /** make index.html as overview */
  WriteHTML(Histos->getPlotsDir() + run_number_, telescopeID);

    getTime(now1, endProg);
    getTime(now2, allProg);
    /** print total events */
    cout << "=======================\n";
    cout << "Total events: " << nEntries << endl;
    cout << "Start       : " << setprecision(2) << fixed << startProg << " seconds\n";
    cout << "Loop        : " << setprecision(2) << fixed << loop << " seconds\n";
    cout << "End         : " << setprecision(2) << fixed << endProg << " seconds\n";
    cout << "All         : " << setprecision(2) << fixed << allProg << " seconds\n";
    cout <<"=======================\n";
}


/** ============================
 AUXILIARY FUNCTIONS
 =================================*/
void PLTAnalysis::SinglePlaneStudies(){

    if ((telescopeID == 1) || (telescopeID == 2)){
        int n_events = TestPlaneEfficiencySilicon(in_file_name_, out_f, run_number_, telescopeID);
        for (uint8_t iplane = 1; iplane != 5; iplane++){
            cout << "Going to call TestPlaneEfficiency " << iplane << endl;
            TestPlaneEfficiency(in_file_name_, out_f, run_number_, iplane, n_events, telescopeID);
        }
    }
}

float PLTAnalysis::getTime(float now, float & time){

    time += (clock() - now) / CLOCKS_PER_SEC;
    return time;
}

void PLTAnalysis::WriteTrackingTree(){

    /** first clear all vectors */
    FW->clearVectors();

    FW->setHitPlaneBits(uint16_t(FR->HitPlaneBits()) );
    FW->setNTracks(uint8_t(FR->NTracks()) );
    FW->setTotalClusters(uint8_t(FR->NClusters()) );

    for (uint8_t iplane = 0; iplane != FR->NPlanes(); ++iplane) {
      PLTPlane * Plane = FR->Plane(iplane);
      FW->setNHits(iplane, uint16_t(Plane->NHits()) );
      FW->setNClusters(iplane, uint8_t(Plane->NClusters()));
      for (size_t icluster = 0; icluster != Plane->NClusters(); icluster++) {
        PLTCluster * Cluster = Plane->Cluster(icluster);
        FW->setClusterPos(iplane, Cluster->SeedHit()->Column(), Cluster->SeedHit()->Row() );
        FW->setClusterPosTel(iplane, Cluster->TX() , Cluster->TY());
        FW->setClusterPosLocal(iplane, Cluster->LX() , Cluster->LY());
        FW->setClusterCharge(iplane, Cluster->Charge());
        FW->setClusterSize(iplane, int(Cluster->NHits()));
      }
    }
    FW->set_dut_tracks(DiaZ); /** set extrapolated position of the track at the diamond position */
    if (FR->NTracks() > 0){
        PLTTrack * Track = FR->Track(0);
        FW->setChi2(Track->Chi2(), Track->Chi2X(), Track->Chi2Y() );
        FW->setAngle(Track->fAngleX, Track->fAngleY);
        for (uint8_t iplane = 0; iplane != FR->NPlanes(); ++iplane) {
            PLTPlane * Plane = FR->Plane(iplane);
            for (size_t icluster = 0; icluster != Plane->NClusters(); icluster++) {
                PLTCluster * Cluster = Plane->Cluster(icluster);
                float xl = FR->GetAlignment()->TtoLX(Track->ExtrapolateX(Cluster->TZ()), Track->ExtrapolateY(Cluster->TZ()), Cluster->Channel(), iplane);
                float yl = FR->GetAlignment()->TtoLY(Track->ExtrapolateX(Cluster->TZ()), Track->ExtrapolateY(Cluster->TZ()), Cluster->Channel(), iplane);
                FW->setResidualXY(iplane, xl - Cluster->LX(), yl - Cluster->LY() );
                FW->setResidual(iplane, float(tel::distance(make_pair(xl, yl), make_pair(Cluster->LX(), Cluster->LY()))) );
                FW->setTrackPos(iplane, Track->ExtrapolateX(Plane->TZ()), Track->ExtrapolateY(Plane->TZ()) );
            }
            FW->setSResidual(iplane, Plane->NClusters() == 1);
        }
    }
    else {
        FW->setChi2(-999, -999, -999);
        FW->setAngle(-999, -999);
    }

    FW->fillTree();
}
void PLTAnalysis::MakeAvgPH(){

    for (uint16_t i = 0; i != Histos->NRoc(); ++i) {
        for (uint16_t j = 0; j != 4; ++j) {
            Histos->AveragePH()[i][j]->Set(NGraphPoints+1);
            Histos->AveragePH()[i][j]->SetPoint(NGraphPoints, ThisTime - TimeWidth/2, Histos->dAveragePH()[i][j]);
            Histos->AveragePH()[i][j]->SetPointError(NGraphPoints, TimeWidth/2, Histos->dAveragePH()[i][j]/sqrt((float) Histos->nAveragePH()[i][j]));
            if (verbose == 1)
                printf("AvgCharge: %i %i N:%9i : %13.3E\n", i, j, Histos->nAveragePH()[i][j], Histos->dAveragePH()[i][j]);
            Histos->dAveragePH()[i][j] = 0;
            Histos->nAveragePH()[i][j] = 0;
        }
    }
    NGraphPoints++;
}
void PLTAnalysis::DrawTracks(){

    static int ieventdraw = 0;
    if (ieventdraw < 20) {
        auto hp = uint16_t(FR->HitPlaneBits());
        if (hp == pow(2, FR->NPlanes() ) -1){
            FR->DrawTracksAndHits(TString::Format(Histos->getOutDir() + "/Tracks_Ev%i.png", ++ieventdraw).Data() );
            if (ieventdraw == 20) cout << endl;
        }
    }
}
void PLTAnalysis::FillPHHistos(uint8_t iplane, PLTCluster * Cluster){

    if (iplane < Histos->NRoc() ) {
        /** fill pulse height histo for all*/
        Histos->PulseHeight()[iplane][0]->Fill(Cluster->Charge());
        Histos->PulseHeightLong()[iplane][0]->Fill(Cluster->Charge());

        /** average pulse heights */
        PLTU::AddToRunningAverage(Histos->dAveragePH2D()[iplane][Cluster->SeedHit()->Column()][ Cluster->SeedHit()->Row()], Histos->nAveragePH2D()[iplane][Cluster->SeedHit()->Column()][ Cluster->SeedHit()->Row()], Cluster->Charge());
        PLTU::AddToRunningAverage(Histos->dAveragePH()[iplane][0], Histos->nAveragePH()[iplane][0], Cluster->Charge());

        /** fill pulse height histo one pix */
        if (Cluster->NHits() == 1) {
            Histos->PulseHeight()[iplane][1]->Fill(Cluster->Charge());
            Histos->PulseHeightLong()[iplane][1]->Fill(Cluster->Charge());
            PLTU::AddToRunningAverage(Histos->dAveragePH()[iplane][1], Histos->nAveragePH()[iplane][1], Cluster->Charge());
        }

        /** fill pulse height histo two pix */
        else if (Cluster->NHits() == 2) {
            Histos->PulseHeight()[iplane][2]->Fill(Cluster->Charge());
            Histos->PulseHeightLong()[iplane][2]->Fill(Cluster->Charge());
            PLTU::AddToRunningAverage(Histos->dAveragePH()[iplane][2], Histos->nAveragePH()[iplane][2], Cluster->Charge());
        }
        /** fill pulse height histo >3 pix */
        else if (Cluster->NHits() >= 3) {
            Histos->PulseHeight()[iplane][3]->Fill(Cluster->Charge());
            Histos->PulseHeightLong()[iplane][3]->Fill(Cluster->Charge());
            PLTU::AddToRunningAverage(Histos->dAveragePH()[iplane][3], Histos->nAveragePH()[iplane][3], Cluster->Charge());
        }
    }
}
void PLTAnalysis::FillOccupancyHiLo(PLTCluster * Cluster){

    /** fill high occupancy */
    if (Cluster->Charge() > 50000)
        for (size_t ihit = 0; ihit != Cluster->NHits(); ++ihit)
            Histos->OccupancyHighPH()[Cluster->ROC()]->Fill( Cluster->Hit(ihit)->Column(), Cluster->Hit(ihit)->Row() );

    /** fill low occupancy */
    else if (Cluster->Charge() > 10000 && Cluster->Charge() < 40000)
        for (size_t ihit = 0; ihit != Cluster->NHits(); ++ihit)
            Histos->OccupancyLowPH()[Cluster->ROC()]->Fill( Cluster->Hit(ihit)->Column(), Cluster->Hit(ihit)->Row() );
}
void PLTAnalysis::FillOccupancy(PLTPlane * Plane){

    for (size_t ihit = 0; ihit != Plane->NHits(); ++ihit) {
        PLTHit * Hit = Plane->Hit(ihit);

        if (Hit->ROC() < Histos->NRoc() ) Histos->Occupancy()[Hit->ROC()]->Fill(Hit->Column(), Hit->Row());
        else cerr << "Oops, ROC >= NROC?" << endl;
    }
}
void PLTAnalysis::FillOfflinePH(PLTTrack * Track, PLTCluster * Cluster){

    if ((fabs(Track->fAngleX) < 0.01) && (fabs(Track->fAngleY) < 0.01)){
						Histos->PulseHeightOffline()[Cluster->ROC()][0]->Fill(Cluster->Charge());

        if (Cluster->NHits() == 1)
            Histos->PulseHeightOffline()[Cluster->ROC()][1]->Fill(Cluster->Charge());
        else if (Cluster->NHits() == 2)
            Histos->PulseHeightOffline()[Cluster->ROC()][2]->Fill(Cluster->Charge());
        else if (Cluster->NHits() >= 3)
            Histos->PulseHeightOffline()[Cluster->ROC()][3]->Fill(Cluster->Charge());
    }
}

vector<float> * PLTAnalysis::getDiaZPositions(){

    auto * tmp = new vector<float>;
    for (uint8_t i_dut(0); i_dut < GetNDUTs(); i_dut++){
      float pos = UseDigitalCalibration() ? FR->GetAlignment()->LZ(1, 4 + i_dut) : tel::Config::dia_z_pos_.at(i_dut);
      cout << "z-position of DUT " << int(i_dut) << ": " << pos << endl;
      tmp->push_back(pos);
    }
  return tmp;
}
