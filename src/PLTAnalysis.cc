#include "PLTAnalysis.h"

using namespace std;

PLTAnalysis::PLTAnalysis(string const inFileName, TFile * Out_f,  TString const runNumber, uint8_t const TelescopeID):
    telescopeID(TelescopeID), InFileName(inFileName), RunNumber(runNumber),
    now1(clock()), now2(clock()), loop(0), startProg(0), endProg(0), allProg(0), averTime(0),
    TimeWidth(20000), StartTime(0), NGraphPoints(0),
    PHThreshold(3e5)
{
    out_f = Out_f;
    /** set up root */
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;
    /** single plane studies */
    SinglePlaneStudies();
    /** init file reader */
    InitFileReader();
    nEntries = ((PSIRootFileReader*) FR)->fTree->GetEntries();
    /** apply masking */
    FR->ReadPixelMask(GetMaskingFilename(telescopeID));
    /** init histos */
    Histos = new RootItems(telescopeID, RunNumber);
    cout << "Output directory: " << Histos->getOutDir() << endl;
    /** init file writer */
    if (telescopeID == 7) FW = new FileWriterTracking(InFileName, telescopeID, FR);
}

PLTAnalysis::~PLTAnalysis()
{
    FW->saveTree();
    delete FR;
    delete FW;
    delete Histos;
}


/** ============================
 EVENT LOOP
 =================================*/
 void PLTAnalysis::EventLoop(){

    getTime(now1, startProg);
    now1 = clock();
    for (uint32_t ievent = 0; FR->GetNextEvent() >= 0; ++ievent) {

        ThisTime = ievent;

        MeasureSpeed(ievent);
        PrintProcess(ievent);

        /** file writer */
        if (GetUseRootInput(telescopeID)) WriteTrackingTree();

        /** fill coincidence map */
        Histos->CoincidenceMap()->Fill(FR->HitPlaneBits() );

        /** make average pulseheight maps*/
        if (ThisTime - (StartTime + NGraphPoints * TimeWidth) > TimeWidth)
            MakeAvgPH();

        /** draw tracks if there is more than one hit*/
        DrawTracks();

        /** loop over the planes */
        for (size_t iplane = 0; iplane != FR->NPlanes(); ++iplane) {

            PLTPlane * Plane = FR->Plane(iplane);

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

		if (telescopeID == 7 && FR->NTracks() == 1 && FR->Track(0)->NClusters() == Histos->NRoc()) {

            do_slope = true;
            for (uint8_t i_rocs(0); i_rocs != Histos->NRoc(); i_rocs++)
                if (FR->Track(0)->Cluster(i_rocs)->Charge() > PHThreshold){
                    do_slope = false;
                    break;
                }

            if (do_slope) {

				PLTTrack * Track = FR->Track(0);

                /** fill chi2 histos */
                Histos->Chi2()->Fill(Track->Chi2());
                Histos->Chi2X()->Fill(Track->Chi2X());
                Histos->Chi2Y()->Fill(Track->Chi2Y());

                /** fill slope histos */
                Histos->TrackSlopeX()->Fill(Track->fSlopeX);
                Histos->TrackSlopeY()->Fill(Track->fSlopeY);

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
	WriteHTML(Histos->getPlotsDir() + RunNumber, GetCalibrationFilename(telescopeID), telescopeID);

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
        int n_events = TestPlaneEfficiencySilicon(InFileName, out_f, RunNumber, telescopeID);
        for (uint8_t iplane = 1; iplane != 5; iplane++){
            cout << "Going to call TestPlaneEfficiency " << iplane << endl;
            TestPlaneEfficiency(InFileName, out_f, RunNumber, iplane, n_events, telescopeID);
        }
    }
}
void PLTAnalysis::InitFileReader(){

    if (GetUseRootInput(telescopeID)){
        FR = new PSIRootFileReader(InFileName, GetCalibrationFilename(telescopeID), GetAlignmentFilename(telescopeID),
            GetNumberOfROCS(telescopeID), GetUseGainInterpolator(telescopeID), GetUseExternalCalibrationFunction(telescopeID));
    }
    else {
        FR = new PSIBinaryFileReader(InFileName, GetCalibrationFilename(telescopeID), GetAlignmentFilename(telescopeID),
            GetNumberOfROCS(telescopeID), GetUseGainInterpolator(telescopeID), GetUseExternalCalibrationFunction(telescopeID));
        ((PSIBinaryFileReader*) FR)->CalculateLevels(10000, Histos->getOutDir());
    }
    FR->GetAlignment()->SetErrors(telescopeID);
    FILE * f = fopen("MyGainCal.dat", "w");
    FR->GetGainCal()->PrintGainCal(f);
    fclose(f);
}
float PLTAnalysis::getTime(float now, float & time){

    time += (clock() - now) / CLOCKS_PER_SEC;
    return time;
}
void PLTAnalysis::PrintProcess(uint32_t ievent){

    if (ievent % 10 == 0 && ievent >= 20000){
        if (ievent != 0) cout << "\e[A\r";
        cout << "Processing event: " << setw(7) << setfill('0') << ievent << endl;
        if (speed) cout << "time left: " << setprecision(2) << fixed << (nEntries - ievent) / speed << "      ";
        else cout << "time left: ???";
    }
}
void PLTAnalysis::MeasureSpeed(uint32_t ievent){

    if (ievent == 10000) now = clock();
    if (ievent % 10000 == 0 && ievent >= 20000){
        speed = (ievent - 10000) / getTime(now, averTime);
        now = clock();
    }
}
void PLTAnalysis::WriteTrackingTree(){

    /** first clear all vectors */
    FW->clearVectors();

    FW->setHitPlaneBits(FR->HitPlaneBits() );
    FW->setNTracks(FR->NTracks() );
    FW->setNClusters(FR->NClusters() );

    if (FR->NTracks() > 0){
        PLTTrack * Track = FR->Track(0);

        FW->setChi2(Track->Chi2() );
        FW->setChi2X(Track->Chi2X() );
        FW->setChi2Y(Track->Chi2Y() );
        FW->setSlopeX(Track->fSlopeX);
        FW->setSlopeY(Track->fSlopeY);
        FW->setDia1TrackX(Track->ExtrapolateX(PLTU::DIA1Z));
        FW->setDia1TrackY(Track->ExtrapolateY(PLTU::DIA1Z));
        FW->setDia2TrackX(Track->ExtrapolateX(PLTU::DIA2Z));
        FW->setDia2TrackY(Track->ExtrapolateY(PLTU::DIA2Z));
        FW->setDistDia1(Track->ExtrapolateX(PLTU::DIA1Z), Track->ExtrapolateY(PLTU::DIA1Z));
        FW->setDistDia2(Track->ExtrapolateX(PLTU::DIA2Z), Track->ExtrapolateY(PLTU::DIA2Z));
    }
    for (size_t iplane = 0; iplane != FR->NPlanes(); ++iplane)
        for (size_t icluster = 0; icluster != FR->Plane(iplane)->NClusters(); ++icluster)
            FW->setChargeAll(iplane, FR->Plane(iplane)->Cluster(icluster)->Charge() );

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
        uint16_t hp = FR->HitPlaneBits();
        if (hp != 0 && hp != 1 && hp != 2 && hp != 4 && hp != 8){
            FR->DrawTracksAndHits(TString::Format(Histos->getOutDir() + "/Tracks_Ev%i.gif", ++ieventdraw).Data() );
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

    if ((fabs(Track->fSlopeX) < 0.01) && (fabs(Track->fSlopeY) < 0.01)){
						Histos->PulseHeightOffline()[Cluster->ROC()][0]->Fill(Cluster->Charge());

        if (Cluster->NHits() == 1)
            Histos->PulseHeightOffline()[Cluster->ROC()][1]->Fill(Cluster->Charge());
        else if (Cluster->NHits() == 2)
            Histos->PulseHeightOffline()[Cluster->ROC()][2]->Fill(Cluster->Charge());
        else if (Cluster->NHits() >= 3)
            Histos->PulseHeightOffline()[Cluster->ROC()][3]->Fill(Cluster->Charge());
    }
}
