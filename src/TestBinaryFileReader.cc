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
    float now1(clock()), now3(clock()), loop(0), startProg(0), endProg(0), allProg(0);
    const uint16_t NROC = GetNumberOfROCS(telescopeID);
    const uint32_t ph_threshold(300000);
    bool do_slope;
    /** set output directory and gStyle */
    TString const PlotsDir = "plots/";
    TString const OutDir = PlotsDir + RunNumber + "/";
    cout << OutDir << endl;
    gStyle->SetOptStat(0);


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

    /** apply masking */
    FR->ReadPixelMask(GetMaskingFilename(telescopeID));


    /** ============================
     Prepare Root Items
     =================================*/
    RootItems RootItems(telescopeID, RunNumber);

    float_t  onepc[NROC];
    float_t  twopc[NROC];
    float_t threepc[NROC];


    /** ============================
     Times for Counting
     ================================= */
    int const TimeWidth = 20000;
    int NGraphPoints = 0;
    int const StartTime = 0;
    int ThisTime;


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

        /** print progress */
        if (ievent % 10000 == 0){
            cout << "Processing event: " << ievent << endl;
            cout << "elapsed time: " << setprecision(2) << fixed << float((clock() - now)) / CLOCKS_PER_SEC << endl;
            now = clock();
        }

        /** file writer */
        if (GetUseRootInput(telescopeID)){

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
            if (hp != 0 && hp != 1 && hp != 2 && hp != 4 && hp != 8)
                FR->DrawTracksAndHits(TString::Format(OutDir + "/Tracks_Ev%i.gif", ++ieventdraw).Data() );
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
                        onepc[iplane]++;
                        RootItems.PulseHeightLong()[iplane][1]->Fill(Cluster->Charge());
                        PLTU::AddToRunningAverage(RootItems.dAveragePH()[iplane][1], RootItems.nAveragePH()[iplane][1], Cluster->Charge());
                    }

                    /** fill pulse height histo two pix */
                    else if (Cluster->NHits() == 2) {
                        RootItems.PulseHeight()[iplane][2]->Fill(Cluster->Charge());
                        twopc[iplane]++;
                        RootItems.PulseHeightLong()[iplane][2]->Fill(Cluster->Charge());
                        PLTU::AddToRunningAverage(RootItems.dAveragePH()[iplane][2], RootItems.nAveragePH()[iplane][2], Cluster->Charge());
                    }
                    /** fill pulse height histo >3 pix */
                    else if (Cluster->NHits() >= 3) {
                        RootItems.PulseHeight()[iplane][3]->Fill(Cluster->Charge());
                        threepc[iplane]++;
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

  TCanvas Can;
  Can.cd();

  out_f->cd();

    RootItems.FormatPHHisto(RootItems.PulseHeight());
    RootItems.FormatPHHisto(RootItems.PulseHeightLong());
    RootItems.FormatPHHisto(RootItems.PulseHeightOffline());
    RootItems.FormatLegendPH();

  for (int iroc = 0; iroc != NROC; ++iroc) {

    /** Draw Occupancy histograms */
    RootItems.Occupancy()[iroc]->SetMinimum(0);
    //RootItems.Occupancy()[iroc].SetAxisRange(12,38,"X");
    //RootItems.Occupancy()[iroc].SetAxisRange(39,80,"Y");
    RootItems.Occupancy()[iroc]->Draw("colz");
    Can.SaveAs( OutDir+TString(RootItems.Occupancy()[iroc]->GetName()) + ".gif");
    RootItems.Occupancy()[iroc]->Write();

    TH1F* hOccupancy1DZ = PLTU::HistFrom2D(RootItems.Occupancy()[iroc]);
    Can.cd();
    hOccupancy1DZ->Draw("hist");
    if (hOccupancy1DZ->GetEntries() > 0) {
      Can.SetLogy(1);
    }
    Can.SaveAs(OutDir+TString(hOccupancy1DZ->GetName()) + ".gif");
    hOccupancy1DZ->Write();
    Can.SetLogy(0);

    // Grab the quantile you're interested in here
    Double_t QProbability[1] = { 0.95 }; // Quantile positions in [0, 1]
    Double_t QValue[1];                  // Quantile values
    hOccupancy1DZ->GetQuantiles(1, QValue, QProbability);
    if(QValue[0] > 1 && RootItems.Occupancy()[iroc]->GetMaximum() > QValue[0]) {
      RootItems.Occupancy()[iroc]->SetMaximum(QValue[0]);
    }
    Can.cd();
    RootItems.Occupancy()[iroc]->Draw("colz");
    Can.SaveAs( OutDir+Form("Occupancy_ROC%i_Quantile.gif", iroc) );
    delete hOccupancy1DZ;

    Can.cd();
    hOccupancy1DZ = PLTU::HistFrom2D(RootItems.Occupancy()[iroc], 0, QValue[0], TString::Format("Occupancy1DZ_ROC%i_Quantile", iroc), 20);
    hOccupancy1DZ->Draw("hist");
    Can.SaveAs(OutDir+TString(hOccupancy1DZ->GetName()) + ".gif");
    delete hOccupancy1DZ;


    // Get 3x3 efficiency hists and draw
    TH2F* h3x3 = PLTU::Get3x3EfficiencyHist(*RootItems.Occupancy()[iroc], 0, 51, 0, 79);
    h3x3->SetTitle( TString::Format("Occupancy Efficiency 3x3 ROC%i", iroc) );
    Can.cd();
    h3x3->SetMinimum(0);
    h3x3->SetMaximum(3);
    h3x3->Draw("colz");
    Can.SaveAs(OutDir+TString(h3x3->GetName()) + ".gif");

    Can.cd();
    TH1F* h3x3_1DZ = PLTU::HistFrom2D(h3x3, "", 50);
    h3x3_1DZ->Draw("hist");
    Can.SaveAs(OutDir+TString(h3x3_1DZ->GetName()) + ".gif");
    delete h3x3;

    /** clusters per event */
    RootItems.DrawSaveTH1F(RootItems.nClusters(), iroc, Can, "Number of clusters per event", "Events");

    /** hits per cluster */
    RootItems.DrawSaveTH1F(RootItems.nHitsPerCluster(), iroc, Can, "Number of hits per cluster", "Number of Clusters");

    Can.cd();
    RootItems.OccupancyHighPH()[iroc]->SetMinimum(0);
    RootItems.OccupancyHighPH()[iroc]->Draw("colz");
    Can.SaveAs( OutDir+TString(RootItems.OccupancyHighPH()[iroc]->GetName()) + ".gif");

    RootItems.OccupancyLowPH()[iroc]->SetMinimum(0);
    RootItems.OccupancyLowPH()[iroc]->Draw("colz");
    Can.SaveAs( OutDir+TString(RootItems.OccupancyLowPH()[iroc]->GetName()) + ".gif");

//    float_t oneovertwo[iroc],oneoverthree[iroc],twooverthree[iroc];  //unused?
//
//    oneovertwo[iroc] = onepc[iroc]/twopc[iroc];
//    oneoverthree[iroc] = onepc[iroc]/threepc[iroc];
//    twooverthree[iroc] = twopc[iroc]/threepc[iroc]; //unused?

    /** ============================
     Draw the Pulse Heights
     ================================= */
    gStyle->SetOptStat(0);

    /** standard */
    RootItems.ClearLegendsPH();
    RootItems.FillLegendsPH(iroc, RootItems.PulseHeight());
    RootItems.DrawSavePH(iroc, RootItems.PulseHeight(), "Pulse Height ROC%i", "PulseHeight_ROC%i.gif");
//    TLegend lRatio(0.75, 0.1, 0.90, 0.4, "Ratio:");
//    lRatio.SetTextAlign(11);
//    lRatio.SetFillStyle(0);
//    lRatio.SetBorderSize(0);
//	  leg_mean.AddEntry( "oneovertwo", TString::Format("%8.0f", oneovertwo[iroc])+" 1pix/2pix");
//    lRatio.AddEntry( "oneoverthree", TString::Format("%8.0f", oneoverthree, "")+" 1 over 3");
//    lRatio.AddEntry( "twooverthree", TString::Format("%8.0f", twooverthree, "")+" 2 over 3");
//    lRatio.Draw("same");
    /** offline */
    RootItems.ClearLegendsPH();
    RootItems.FillLegendsPH(iroc, RootItems.PulseHeightOffline());
    RootItems.DrawSavePH(iroc, RootItems.PulseHeightOffline(), "Pulse Height Offline ROC%i", "PulseHeightOffline_ROC%i.gif");

    /** long */
    RootItems.ClearLegendsPH();
    RootItems.FillLegendsPH(iroc, RootItems.PulseHeightLong());
    RootItems.DrawSavePH(iroc, RootItems.PulseHeightLong(), "Pulse Height Long ROC%i", "PulseHeightLong_ROC%i.gif");

    /** average pulse height */
    RootItems.DrawSaveAvPH(iroc);

    // Use AvgPH2D to draw PH 2D maps
    TString Name = TString::Format("PulseHeightAvg2D_ROC%i", iroc);
    TH2F hPulseHeightAvg2D(Name, Name, PLTU::NCOL, PLTU::FIRSTCOL, PLTU::LASTCOL, PLTU::NROW, PLTU::FIRSTROW, PLTU::LASTROW);
    for (int icol = 0; icol != PLTU::NCOL; ++icol) {
      for (int irow = 0; irow != PLTU::NROW; ++irow) {
        hPulseHeightAvg2D.SetBinContent(icol+1, irow+1, RootItems.dAveragePH2D()[iroc][icol][irow]);
//        hPulseHeightAvg2D.SetBinContent(icol+1, irow+1, AvgPH2D[iroc][icol][irow]);
      }
    }
    Can.cd();
    hPulseHeightAvg2D.SetMinimum(0);
    hPulseHeightAvg2D.SetMaximum(100000);
    hPulseHeightAvg2D.Draw("colz");
    Can.SaveAs(OutDir+hPulseHeightAvg2D.GetName() + ".gif");


    /** residuals */
    RootItems.DrawSaveResidual(iroc, RootItems.Residual());
    RootItems.DrawSaveResidual(iroc, RootItems.ResidualXdY());
    RootItems.DrawSaveResidual(iroc, RootItems.ResidualYdX());
    RootItems.DrawSaveResidualProj(iroc, RootItems.Residual(), "X");
    RootItems.DrawSaveResidualProj(iroc, RootItems.Residual(), "Y");

  } // end of loop over ROCs

    /** draw and save coincidence map */
    RootItems.PrepCoincidenceHisto();
    RootItems.DrawSaveCoincidence();

    /** draw tracking slopes and chi2 */
    gStyle->SetOptFit(101);
    Can.cd();
    RootItems.FitSlope(RootItems.TrackSlopeX() );
    RootItems.TrackSlopeX()->Draw();
    RootItems.LegendSlope(RootItems.TrackSlopeX() );
    Can.SaveAs(OutDir + "TrackSlopeX.gif");
    RootItems.TrackSlopeX()->Write();

    Can.cd();
    RootItems.FitSlope(RootItems.TrackSlopeY() );
    RootItems.TrackSlopeY()->Draw();
    RootItems.LegendSlope(RootItems.TrackSlopeY() );
    Can.SaveAs(OutDir+"TrackSlopeY.gif");
    RootItems.TrackSlopeY()->Write();

    RootItems.DrawSaveChi2(RootItems.Chi2(), "Chi2");
    RootItems.DrawSaveChi2(RootItems.Chi2X(), "Chi2X");
    RootItems.DrawSaveChi2(RootItems.Chi2Y(), "Chi2Y");

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
