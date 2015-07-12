#include "TestBinaryFileReader.h"
#include "RootItems.h"
#include "FileWriterTracking.h"

#define verbose 0

using namespace std;

/** ============================
 RUN DEFAULT ANALYSIS
 =================================*/

float getTime(float now, float & time){
    time += (clock() - now) / CLOCKS_PER_SEC;
    return time;
}

int TestPSIBinaryFileReader (string const InFileName, TFile * out_f,  TString const RunNumber,
                             int telescopeID)
{

    /** ============================
     Constants
     =================================*/
    /** measure elapsed time */
    float now1(clock()), now3(clock()), loop(0), startProg(0), endProg(0), allProg(0), averTime(0), speed;
    /** times for counting */
    int const TimeWidth = 20000;
    int NGraphPoints = 0;
    int const StartTime = 0;
    int ThisTime;
    /** miscellaneous */
    const uint16_t NROC = GetNumberOfROCS(telescopeID);
    const uint32_t ph_threshold(300000);
    bool do_slope;
    /** set output directory and gStyle */
    TString const PlotsDir = "plots/";
    TString const OutDir = PlotsDir + RunNumber + "/";
    cout << OutDir << endl;
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;


    /** ============================
     Single Plane Studies
     =================================*/
    if ((telescopeID == 1) || (telescopeID == 2)){
        int n_events = TestPlaneEfficiencySilicon(InFileName, out_f, RunNumber, telescopeID);
        for (int iplane = 1; iplane != 5; iplane++){
            cout << "Going to call TestPlaneEfficiency " << iplane << endl;
            TestPlaneEfficiency(InFileName, out_f, RunNumber, iplane, n_events, telescopeID);
        }
    }


    /** ============================
     Initialize File Reader
     =================================*/
    PSIFileReader * FR;

    if (GetUseRootInput(telescopeID)){
        FR = new PSIRootFileReader(InFileName, GetCalibrationFilename(telescopeID), GetAlignmentFilename(telescopeID),
            GetNumberOfROCS(telescopeID), GetUseGainInterpolator(telescopeID), GetUseExternalCalibrationFunction(telescopeID));
    }
    else {
        FR = new PSIBinaryFileReader(InFileName, GetCalibrationFilename(telescopeID), GetAlignmentFilename(telescopeID),
            GetNumberOfROCS(telescopeID), GetUseGainInterpolator(telescopeID), GetUseExternalCalibrationFunction(telescopeID));
        ((PSIBinaryFileReader*) FR)->CalculateLevels(10000, OutDir);
    }
    FR->GetAlignment()->SetErrors(telescopeID);
    FILE * f = fopen("MyGainCal.dat", "w");
    FR->GetGainCal()->PrintGainCal(f);
    fclose(f);
    uint32_t nEntries = ((PSIRootFileReader*) FR)->fTree->GetEntries();

    /** apply masking */
    FR->ReadPixelMask(GetMaskingFilename(telescopeID));


    /** ============================
     Prepare Root Items
     =================================*/
    RootItems RootItems(telescopeID, RunNumber);


    /** ============================
     Initialize File Writer
     =================================*/
    FileWriterTracking * FW;
    if (GetUseRootInput(telescopeID) && (telescopeID==7))
        FW = new FileWriterTracking(InFileName, telescopeID, FR);


    getTime(now1, startProg);
    now1 = clock();
    /** ============================
     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
     EVENT LOOP
     xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
     ================================= */
    uint64_t indexi = -1;
    float now = clock();
    for (uint32_t ievent = 0; FR->GetNextEvent() >= 0; ++ievent) {

        indexi++;
        ThisTime = ievent;

        /** print process */
        if (ievent % 10 == 0 && ievent >= 20000){
            if (ievent != 0) cout << "\e[A\r";
            cout << "Processing event: " << setw(7) << setfill('0') << ievent << endl;
            if (speed) cout << "time left: " << setprecision(2) << fixed << (nEntries - ievent) / speed << "      ";
            else cout << "time left: ???";
        }

        /** measure speed */
        if (ievent % 10000 == 0 && ievent >= 20000){
            speed = ievent / getTime(now, averTime);
            now = clock();
        }

        /** file writer */
        if (GetUseRootInput(telescopeID)){

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

        /** fill coincidence map */
        RootItems.CoincidenceMap()->Fill(FR->HitPlaneBits());

        /** make average pulseheight maps*/
        if (ThisTime - (StartTime + NGraphPoints * TimeWidth) > TimeWidth) {
            for (uint16_t i = 0; i != NROC; ++i) {
                for (uint16_t j = 0; j != 4; ++j) {
                    RootItems.AveragePH()[i][j]->Set(NGraphPoints+1);
                    RootItems.AveragePH()[i][j]->SetPoint(NGraphPoints, ThisTime - TimeWidth/2, RootItems.dAveragePH()[i][j]);
                    RootItems.AveragePH()[i][j]->SetPointError(NGraphPoints, TimeWidth/2, RootItems.dAveragePH()[i][j]/sqrt((float) RootItems.nAveragePH()[i][j]));
                    if (verbose == 1)
                        printf("AvgCharge: %i %i N:%9i : %13.3E\n", i, j, RootItems.nAveragePH()[i][j], RootItems.dAveragePH()[i][j]);
                    RootItems.dAveragePH()[i][j] = 0;
                    RootItems.nAveragePH()[i][j] = 0;
                }
            }
          ++NGraphPoints;
        }

        /** draw tracks if there is more than one hit*/
        static int ieventdraw = 0;
        if (ieventdraw < 20) {
            uint16_t hp = FR->HitPlaneBits();
            if (hp != 0 && hp != 1 && hp != 2 && hp != 4 && hp != 8){
                FR->DrawTracksAndHits(TString::Format(OutDir + "/Tracks_Ev%i.gif", ++ieventdraw).Data() );
                if (ieventdraw == 20) cout << endl;
            }
        }

        /** loop over the planes */
        for (size_t iplane = 0; iplane != FR->NPlanes(); ++iplane) {

            PLTPlane* Plane = FR->Plane(iplane);

            /** fill cluster histo */
            RootItems.nClusters()[Plane->ROC()]->Fill(Plane->NClusters());

            for (size_t icluster = 0; icluster != Plane->NClusters(); ++icluster) {
                PLTCluster* Cluster = Plane->Cluster(icluster);

                if (iplane < NROC) {
                    /** fill pulse height histo for all*/
                    RootItems.PulseHeight()[iplane][0]->Fill(Cluster->Charge());
                    RootItems.PulseHeightLong()[iplane][0]->Fill(Cluster->Charge());

                    /** ignore charges above a certain threshold */
                    if (Cluster->Charge() > ph_threshold) continue;

                    /** average pulse heights */
                    PLTU::AddToRunningAverage(RootItems.dAveragePH2D()[iplane][Cluster->SeedHit()->Column()][ Cluster->SeedHit()->Row()], RootItems.nAveragePH2D()[iplane][Cluster->SeedHit()->Column()][ Cluster->SeedHit()->Row()], Cluster->Charge());
                    PLTU::AddToRunningAverage(RootItems.dAveragePH()[iplane][0], RootItems.nAveragePH()[iplane][0], Cluster->Charge());

                    /** fill pulse height histo one pix */
                    if (Cluster->NHits() == 1) {
                        RootItems.PulseHeight()[iplane][1]->Fill(Cluster->Charge());
                        RootItems.PulseHeightLong()[iplane][1]->Fill(Cluster->Charge());
                        PLTU::AddToRunningAverage(RootItems.dAveragePH()[iplane][1], RootItems.nAveragePH()[iplane][1], Cluster->Charge());
                    }

                    /** fill pulse height histo two pix */
                    else if (Cluster->NHits() == 2) {
                        RootItems.PulseHeight()[iplane][2]->Fill(Cluster->Charge());
                        RootItems.PulseHeightLong()[iplane][2]->Fill(Cluster->Charge());
                        PLTU::AddToRunningAverage(RootItems.dAveragePH()[iplane][2], RootItems.nAveragePH()[iplane][2], Cluster->Charge());
                    }
                    /** fill pulse height histo >3 pix */
                    else if (Cluster->NHits() >= 3) {
                        RootItems.PulseHeight()[iplane][3]->Fill(Cluster->Charge());
                        RootItems.PulseHeightLong()[iplane][3]->Fill(Cluster->Charge());
                        PLTU::AddToRunningAverage(RootItems.dAveragePH()[iplane][3], RootItems.nAveragePH()[iplane][3], Cluster->Charge());
                    }
                }


                /** fill hits per cluster histo */
                RootItems.nHitsPerCluster()[Cluster->ROC()]->Fill(Cluster->NHits());

                /** fill high occupancy high */
                if (Cluster->Charge() > 50000)
                    for (size_t ihit = 0; ihit != Cluster->NHits(); ++ihit)
                        RootItems.OccupancyHighPH()[Cluster->ROC()]->Fill( Cluster->Hit(ihit)->Column(), Cluster->Hit(ihit)->Row() );

                /** fill low occupancy */
                else if (Cluster->Charge() > 10000 && Cluster->Charge() < 40000)
                    for (size_t ihit = 0; ihit != Cluster->NHits(); ++ihit)
                        RootItems.OccupancyLowPH()[Cluster->ROC()]->Fill( Cluster->Hit(ihit)->Column(), Cluster->Hit(ihit)->Row() );
            }
            /** fill occupancy histo */
            for (size_t ihit = 0; ihit != Plane->NHits(); ++ihit) {
                PLTHit* Hit = Plane->Hit(ihit);

                if (Hit->ROC() < NROC) RootItems.Occupancy()[Hit->ROC()]->Fill(Hit->Column(), Hit->Row());
                else cerr << "Oops, ROC >= NROC?" << endl;
            }
        }

		if (telescopeID == 7 && FR->NTracks() == 1 && FR->Track(0)->NClusters() == NROC) {

            do_slope = true;
            for (uint8_t i_rocs(0); i_rocs != NROC; i_rocs++)
                if (FR->Track(0)->Cluster(i_rocs)->Charge() > ph_threshold){
                    do_slope = false;
                    break;
                }

            if (do_slope) {

				PLTTrack * Track = FR->Track(0);

                /** fill chi2 histos */
                RootItems.Chi2()->Fill(Track->Chi2());
                RootItems.Chi2X()->Fill(Track->Chi2X());
                RootItems.Chi2Y()->Fill(Track->Chi2Y());

                /** fill slope histos */
                RootItems.TrackSlopeX()->Fill(Track->fSlopeX);
                RootItems.TrackSlopeY()->Fill(Track->fSlopeY);

				/** loop over the clusters */
				for (size_t icluster = 0; icluster < Track->NClusters(); icluster++){

					/** get the ROC in of the cluster and fill the corresponding residual */
					uint8_t ROC = Track->Cluster(icluster)->ROC();
					PLTCluster * Cluster = Track->Cluster(icluster);

                    /** fill residuals */
					RootItems.Residual()[ROC]->Fill(Track->LResidualX(ROC), Track->LResidualY(ROC)); // dX vs dY
					RootItems.ResidualXdY()[ROC]->Fill(Cluster->LX(), Track->LResidualY(ROC));// X vs dY
					RootItems.ResidualYdX()[ROC]->Fill(Cluster->LY(), Track->LResidualX(ROC)); // Y vs dX

                    /** ignore events above a certain threshold */
					if (Cluster->Charge() > ph_threshold) continue;
					//printf("High Charge: %13.3E\n", Cluster->Charge());

					/** fill the offline pulse heights (Track6+|Slope| < 0.01 in x and y ) */
					if ((fabs(Track->fSlopeX) < 0.01) && (fabs(Track->fSlopeY) < 0.01)){
						RootItems.PulseHeightOffline()[Cluster->ROC()][0]->Fill(Cluster->Charge());

						if (Cluster->NHits() == 1)
							RootItems.PulseHeightOffline()[Cluster->ROC()][1]->Fill(Cluster->Charge());
						else if (Cluster->NHits() == 2)
							RootItems.PulseHeightOffline()[Cluster->ROC()][2]->Fill(Cluster->Charge());
						else if (Cluster->NHits() >= 3)
							RootItems.PulseHeightOffline()[Cluster->ROC()][3]->Fill(Cluster->Charge());
					}
				}
            }
		}

	} /** END OF EVENT LOOP */
	cout << endl;
    getTime(now1, loop);
    now1 = clock();

    if (GetUseRootInput(telescopeID)) FW->saveTree();

    delete FR;
    delete FW;

    /** add the last point to the average pulse height graph */
    for (int i = 0; i != NROC; ++i) {
        for (int j = 0; j != 4; ++j) {
            RootItems.AveragePH()[i][j]->Set(NGraphPoints+1);
            RootItems.AveragePH()[i][j]->SetPoint(NGraphPoints, NGraphPoints*TimeWidth + TimeWidth/2, RootItems.dAveragePH()[i][j]);
            RootItems.AveragePH()[i][j]->SetPointError(NGraphPoints, TimeWidth/2, RootItems.dAveragePH()[i][j]/sqrt((float) RootItems.nAveragePH()[i][j]));
            if (verbose == 1)
                printf("AvgCharge: %i %i N:%9i : %13.3E\n", i, j, RootItems.nAveragePH()[i][j], RootItems.dAveragePH()[i][j]);
        }
    }

    out_f->cd();

    RootItems.SaveAllHistos();

    /** make index.html as overview */
	WriteHTML(PlotsDir + RunNumber, GetCalibrationFilename(telescopeID), telescopeID);

    getTime(now1, endProg);
    getTime(now3, allProg);
    /** print total events */
    cout << "=======================\n";
    cout << "Total events: " << indexi + 1 << endl;
    cout << "Start       : " << setprecision(2) << fixed << startProg << " seconds\n";
    cout << "Loop        : " << setprecision(2) << fixed << loop << " seconds\n";
    cout << "End         : " << setprecision(2) << fixed << endProg << " seconds\n";
    cout << "All         : " << setprecision(2) << fixed << allProg << " seconds\n";
    cout <<"=======================\n";

	return 0;
}
