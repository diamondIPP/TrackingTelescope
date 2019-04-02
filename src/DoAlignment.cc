#include <utility>

#include <Utils.h>
#include <DoAlignment.h>
#include "PSIBinaryFileReader.h"
#include "PSIRootFileReader.h"
#include "GetNames.h"
#include "TestPlaneEfficiencySilicon.h"
#include "PLTPlane.h"
#include "PLTAlignment.h"

using namespace std;

Alignment::Alignment(string in_file_name, TString run_number, short telescope_ID):
  TelescopeID(telescope_ID),
  NPlanes(GetNumberOfROCS(telescope_ID)),
  InFileName(std::move(in_file_name)),
  PlotsDir("plots/"),
  OutDir(PlotsDir + run_number),
  AngleThreshold(.01),
  TotResThreshold(.01),
  Now(clock()) {

  gROOT->ProcessLine("gErrorIgnoreLevel = kError;");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(53);

  XAlign.resize(NPlanes, 0);
  YAlign.resize(NPlanes, 0);
  ZAlign.resize(NPlanes, 0);
  RAlign.resize(NPlanes, 0);

  FR = InitFileReader();
  MaxEventNumber = unsigned(dynamic_cast<PSIRootFileReader*>(FR)->fTree->GetEntries());
  ProgressBar = new tel::ProgressBar(MaxEventNumber - 1);
  /** Apply Masking */
  FR->ReadPixelMask(GetMaskingFilename(TelescopeID));
  InitHistograms();
    
  PreAlign();
}

PSIFileReader * Alignment::InitFileReader() {

  PSIFileReader * tmp;
  if (GetUseRootInput(TelescopeID)){
    tmp = new PSIRootFileReader(InFileName, GetCalibrationFilename(TelescopeID), GetAlignmentFilename(TelescopeID), NPlanes, GetUseGainInterpolator(TelescopeID),
      GetUseExternalCalibrationFunction(TelescopeID), false, uint8_t(TelescopeID));
  }
  else {
    tmp = new PSIBinaryFileReader(InFileName, GetCalibrationFilename(TelescopeID), GetAlignmentFilename(TelescopeID), NPlanes, GetUseGainInterpolator(TelescopeID),
      GetUseExternalCalibrationFunction(TelescopeID));
  }
  tmp->GetAlignment()->SetErrors(TelescopeID);
  FILE * f = fopen("MyGainCal.dat", "w");
  tmp->GetGainCal()->PrintGainCal(f);
  fclose(f);
  return tmp;
}

void Alignment::PreAlign() {

  /** coarsely move the two inner planes to minimise the residuals without rotations */
  cout << "=== STARTING PRE-ALIGNMENT ===" << endl;
  FR->ResetFile();
  FR->SetPlanesUnderTest(1, 2); /** ignore inner planes for tracking */;
  ClearHistograms(); /** Reset residual histograms */
  Now = clock();
        
  /** EVENT LOOP */
  for (uint32_t i_event = 0; FR->GetNextEvent() >= 0; ++i_event) {

    cout << i_event << endl;
    if (i_event > MaxEventNumber) break;
    ProgressBar->update(i_event); /** print progress */
    if (FR->NTracks() != 1) continue; /** proceed only if we have exactly one track */

    PLTTrack * Track = FR->Track(0);
    for (unsigned i_plane(1); i_plane < NPlanes - 1; ++i_plane) {

      if (FR->Plane(i_plane)->NClusters() != 1) continue; /** proceed only if there is exactly one cluster */

      PLTCluster * Cluster = FR->Plane(i_plane)->Cluster(0);
      pair<float, float> l_res = Track->GetResiduals(*Cluster, *FR->GetAlignment());
      cout << "Track Residuals for " << i_plane << " (x/y): " << l_res.first << "/" << l_res.second << endl;
      if (fabs(l_res.first) >= 2 or fabs(l_res.second) >= 2) continue; /** only proceed if the residual is smaller than 2mm in x or y */

      hResidual[i_plane].Fill(l_res.first, l_res.second); // dX vs dY
      hResidualXdY[i_plane].Fill(Cluster->LX(), l_res.second); // X vs dY
      hResidualYdX[i_plane].Fill(Cluster->LY(), l_res.first); // Y vs dX
    }
  } // end event loop
  cout << "\nLoop duration:" << (clock() - Now) / CLOCKS_PER_SEC << endl;



  for (unsigned i_plane(1); i_plane < NPlanes - 1; i_plane++){

    XAlign[i_plane] += hResidual[i_plane].GetMean(1);
    YAlign[i_plane] += hResidual[i_plane].GetMean(2);

    cout << "RESIDUALS Plane " << i_plane << ":\t" << setprecision(7) << XAlign[i_plane] << " " << YAlign[i_plane] << endl;
    cout << "RESIDUALS RMSPlane " << i_plane << ":\t" << setprecision(7) << hResidual[i_plane].GetRMS(1) << " " << hResidual[i_plane].GetRMS(2) <<endl;

    /** update the alignment */
    cout << "Before: " << FR->GetAlignment()->LX(1, i_plane) << endl;
    FR->GetAlignment()->AddToLX(1, i_plane, XAlign[i_plane] );
    FR->GetAlignment()->AddToLY(1, i_plane, YAlign[i_plane] );
    cout << "After:  " << FR->GetAlignment()->LX(1, i_plane) << endl;

    SaveResiduals(i_plane);
  }
}

void Alignment::InitHistograms() {

  for (uint8_t i_plane(0); i_plane != NPlanes; ++i_plane){
    hResidual.emplace_back(        Form("Residual_ROC%i", i_plane),    Form("Residual_ROC%i",    i_plane), 200, -.5, .5, 200, -.5, .5);
    hResidualXdY.emplace_back(TH2F(Form("ResidualXdY_ROC%i", i_plane), Form("ResidualXdY_ROC%i", i_plane), 133, -1, 0.995, 100, -.5, .5));
    hResidualYdX.emplace_back(TH2F(Form("ResidualYdX_ROC%i", i_plane), Form("ResidualYdX_ROC%i", i_plane), 201, -1, 1, 100, -.5, .5));
  }
}

void Alignment::ClearHistograms() {

  for (uint8_t i_plane(0); i_plane != NPlanes; ++i_plane){
    hResidual.at(i_plane).Clear();
    hResidualXdY.at(i_plane).Clear();
    hResidualYdX.at(i_plane).Clear();
  }
}

void Alignment::SaveResiduals(unsigned i_plane) {

  TCanvas Can;
  Can.cd();
  // 2D Residuals
  hResidual[i_plane].SetContour(1024);
  hResidual[i_plane].Draw("colz");
  Can.SaveAs( OutDir + "/" + TString(hResidual[i_plane].GetName()) + ".png");
  // Residual X-Projection
  gStyle->SetOptStat(1111);
  hResidual[i_plane].ProjectionX()->Draw();
  Can.SaveAs( OutDir + "/" + TString(hResidual[i_plane].GetName()) + "_X.png");
  // Residual Y-Projection
  hResidual[i_plane].ProjectionY()->Draw();
  Can.SaveAs( OutDir + "/" + TString(hResidual[i_plane].GetName()) + "_Y.png");
  // 2D Residuals X/dY
  hResidualXdY[i_plane].SetContour(1024);
  hResidualXdY[i_plane].Draw("colz");
  Can.SaveAs( OutDir + "/" + TString(hResidualXdY[i_plane].GetName()) + ".png");
  // 2D Residuals Y/dX
  hResidualYdX[i_plane].SetContour(1024);
  hResidualYdX[i_plane].Draw("colz");
  Can.SaveAs( OutDir + "/" + TString(hResidualYdX[i_plane].GetName()) + ".png");
}

/** DoAlignment: Produce alignment constants and save them to NewAlignment.dat */
int DoAlignment (string const InFileName, TString const RunNumber, int telescopeID) {
    gROOT->ProcessLine("gErrorIgnoreLevel = kError;");
    TString const PlotsDir = "plots/";
    TString const OutDir = PlotsDir + RunNumber;
    float const angle_threshold = 0.01;
    float const tot_res_threshold = 0.01;

    gStyle->SetOptStat(0);
    gStyle->SetPalette(53);

    std::vector<float> x_align;
    std::vector<float> y_align;
    std::vector<float> z_align;
    std::vector<float> r_align;

    for (int i=0; i!=GetNumberOfROCS(telescopeID);i++){
        x_align.push_back(0);
        y_align.push_back(0);
        z_align.push_back(0);
        r_align.push_back(0);
    }

    /** Initialize File Reader*/
    PSIFileReader * FR;

    if (GetUseRootInput(telescopeID)){
        FR = new PSIRootFileReader(InFileName, GetCalibrationFilename(telescopeID), GetAlignmentFilename(telescopeID),
                                   GetNumberOfROCS(telescopeID), GetUseGainInterpolator(telescopeID), GetUseExternalCalibrationFunction(telescopeID), true, telescopeID); // DA: Added telescopeID
    }
    else{
        FR = new PSIBinaryFileReader(InFileName, GetCalibrationFilename(telescopeID), GetAlignmentFilename(telescopeID, true),
                                     GetNumberOfROCS(telescopeID), GetUseGainInterpolator(telescopeID), GetUseExternalCalibrationFunction(telescopeID) );
        ((PSIBinaryFileReader*) FR)->CalculateLevels(OutDir);
    }

    uint32_t stopAt = ((PSIRootFileReader*) FR)->fTree->GetEntries();
//    uint32_t stopAt = 5e4;
    auto * PBar = new tel::ProgressBar(stopAt - 1);

    /** Apply Masking */
    FR->ReadPixelMask(GetMaskingFilename(telescopeID));

    for (int ialign=0; ialign!=2; ialign++){

        for (uint8_t iroc=1; iroc!=GetNumberOfROCS(telescopeID); iroc++){// DA: changed
            FR->GetAlignment()->AddToLX( 1, iroc, x_align[iroc] );
            FR->GetAlignment()->AddToLY( 1, iroc, y_align[iroc] );
            //FR->GetAlignment()->AddToLR( 1, iroc, r_align[iroc] );
        }

        for (int iroc_align = 1; iroc_align != GetNumberOfROCS(telescopeID); ++iroc_align) {

            std::cout << "GOING TO ALIGN: " << iroc_align << std::endl;

            FR->ResetFile();
            FR->SetPlaneUnderTest(iroc_align);// ignore plane for tracking

            /** Prepare Residual histograms
                hResidual:    x=dX / y=dY
                hResidualXdY: x=X  / y=dY
                hResidualYdX: x=Y  / y=dX */
            std::vector<TH2F> hResidual;
            std::vector<TH2F> hResidualXdY;
            std::vector<TH2F> hResidualYdX;

            /** Reset residual histograms */
            hResidual.clear();
            hResidualXdY.clear();
            hResidualYdX.clear();

            for (uint8_t iroc = 0; iroc != GetNumberOfROCS(telescopeID); ++iroc){
                hResidual.push_back(    TH2F(  Form("Residual_ROC%i",iroc),     Form("Residual_ROC%i",iroc), 200, -.5, .5, 200, -.5, .5));
                hResidualXdY.push_back( TH2F(  Form("ResidualXdY_ROC%i",iroc),  Form("ResidualXdY_ROC%i",iroc), 133, -1, 0.995, 100, -.5, .5));
                hResidualYdX.push_back( TH2F(  Form("ResidualYdX_ROC%i",iroc),  Form("ResidualYdX_ROC%i",iroc), 201, -1, 1, 100, -.5, .5));
            }

            float now = clock();
            /** EVENT LOOP */
            for (uint32_t ievent = 0; FR->GetNextEvent() >= 0; ++ievent) {

                cout << ievent << endl;
                if (ievent > stopAt)
                    break;

                /** print progress */
                PBar->update(ievent);

                /** continue only if we have exactly one track */
                if ( !(FR->NTracks()==1) ) continue;

                PLTTrack * Track = FR->Track(0);

                /** continue only if we there is exactly one cluster */
                if ( !(FR->Plane(iroc_align)->NClusters()==1) ) continue;

                float max_charge = -1;
                float h_LX = -9999;
                float h_LY = -9999;

                /** loop over hits in the cluster of the plane under test. Find highest charge hit for residuals*/
                for (uint8_t i=0; i != FR->Plane(iroc_align)->Cluster(0)->NHits(); ++i){

                    PLTHit * Hit = FR->Plane(iroc_align)->Cluster(0)->Hit(i);

                    if (Hit->Charge() > max_charge){
                        max_charge = Hit->Charge();
                        h_LX = Hit->LX();
                        h_LY = Hit->LY();
                    }
                }

                float track_TX = Track->TX(iroc_align); //MR: should be z-position not roc number...
                cout << Track->NClusters() << " " << Track->Cluster(0)->LX() <<", "<< h_LX << " " << h_LY <<", " << track_TX << " " << Track->ExtrapolateX(FR->Plane(iroc_align)->TZ()) << endl;
                float track_TY = Track->TY(iroc_align);

                float track_LX = FR->GetAlignment()->TtoLX( track_TX, track_TY, 1, iroc_align);// Local position of the track in the plane under test
                float track_LY = FR->GetAlignment()->TtoLY( track_TX, track_TY, 1, iroc_align);

                float d_LX =  (track_LX - h_LX);// residuals of track local position and highes charge hit local position
                float d_LY =  (track_LY - h_LY);

                std::cout << "Track LX/LY" << track_LX << " " << track_LY << std::endl;
                // if the residual of track with biggest charge hit is greater than 2mm in x or y, dont take it into account
                if ( !(fabs(d_LX)<2000) || !(fabs(d_LY)<2000) ) continue;// DA: before was <2

                hResidual[iroc_align].Fill( d_LX, d_LY); // dX vs dY

                hResidualXdY[iroc_align].Fill( h_LX, d_LY); // X vs dY

                hResidualYdX[iroc_align].Fill( h_LY, d_LX); // Y vs dX

            } // end event loop
            cout << "\nLoop duration:" << (clock() - now) / CLOCKS_PER_SEC << endl;

            std::cout << "RESIDUALS:\t" << setprecision(7) << hResidual[iroc_align].GetMean(1) << " " << hResidual[iroc_align].GetMean(2) << std::endl;
            std::cout << "RESIDUALS RMS:\t" << setprecision(7) << hResidual[iroc_align].GetRMS(1) << " " << hResidual[iroc_align].GetRMS(2) <<std::endl;

            std::cout << "Before: " << FR->GetAlignment()->LX(1,iroc_align) << std::endl;
            x_align[iroc_align] +=  hResidual[iroc_align].GetMean(1);
            y_align[iroc_align] +=  hResidual[iroc_align].GetMean(2);
            r_align[iroc_align] +=  hResidualXdY[iroc_align].GetCorrelationFactor();
            std::cout << "After:  " << FR->GetAlignment()->LX(1,iroc_align) << std::endl;


            TCanvas Can;
            Can.cd();

            // 2D Residuals
            hResidual[iroc_align].SetContour(1024);
            hResidual[iroc_align].Draw("colz");
            Can.SaveAs( OutDir+"/"+TString(hResidual[iroc_align].GetName()) + ".png");

            // Residual X-Projection
            gStyle->SetOptStat(1111);
            hResidual[iroc_align].ProjectionX()->Draw();
            Can.SaveAs( OutDir+"/"+TString(hResidual[iroc_align].GetName()) + "_X.png");

            // Residual Y-Projection
            hResidual[iroc_align].ProjectionY()->Draw();
            Can.SaveAs( OutDir+"/"+TString(hResidual[iroc_align].GetName()) + "_Y.png");

            // 2D Residuals X/dY
            hResidualXdY[iroc_align].SetContour(1024);
            hResidualXdY[iroc_align].Draw("colz");
            Can.SaveAs( OutDir+"/"+TString(hResidualXdY[iroc_align].GetName()) + ".png");

            // 2D Residuals Y/dX
            hResidualYdX[iroc_align].SetContour(1024);
            hResidualYdX[iroc_align].Draw("colz");
            Can.SaveAs( OutDir+"/"+TString(hResidualYdX[iroc_align].GetName()) + ".png");

//            for (uint8_t i = 0; i != GetNumberOfROCS(telescopeID); i++){
//                std::cout << int(i) << setprecision(7) << " " << x_align[i] << " " << y_align[i] << " " << z_align[i] << " " << r_align[i] <<std::endl;
//            }
        } // end loop over rocs

//        FR->GetAlignment()->WriteAlignmentFile("NewAlignment.dat", FR->NMAXROCS);

    } // end alignment loop




    std::cout << "\n************\nPART TWO!!!!!\n************\n" << std::endl;


    for (int ialign=1; ialign!=15;ialign++){// DA: TODO: use also threshold condition

        std::cout << "BEGIN ITERATION " << ialign << " OUT OF 14" << std::endl;

        FR->ResetFile();
        FR->SetAllPlanes();

        // Prepare Residual histograms
        // hResidual:    x=dX / y=dY
        // hResidualXdY: x=X  / y=dY
        // hResidualYdX: x=Y  / y=dX
        std::vector< TH2F > hResidual;
        std::vector< TH2F > hResidualXdY;
        std::vector< TH2F > hResidualYdX;
        std::vector< TGraph > gResidualXdY;
        std::vector< TGraph > gResidualYdX;

        // Reset residual histograms
        hResidual.clear();
        hResidualXdY.clear();
        hResidualYdX.clear();
        for (int iroc = 0; iroc != GetNumberOfROCS(telescopeID); ++iroc){
            hResidual.push_back( TH2F(  Form("Residual_ROC%i",iroc),
                                        Form("Residual_ROC%i",iroc), 400, -.8, .8, 400, -.8, .8));
            hResidualXdY.push_back( TH2F(  Form("ResidualXdY_ROC%i",iroc),
                                           Form("ResidualXdY_ROC%i",iroc), 35, -0.2, 0.2, 100, -.2, .2));
            hResidualYdX.push_back( TH2F(  Form("ResidualYdX_ROC%i",iroc),
                                           Form("ResidualYdX_ROC%i",iroc), 41, -.2, .2, 100, -.2, .2));
            gResidualXdY.push_back( TGraph() );
            gResidualYdX.push_back( TGraph() );
        }

        // Event Loop
        for (uint32_t ievent = 0; FR->GetNextEvent() >= 0; ++ievent) {

            if (ievent > stopAt)
                break;

            /** print progress */
            PBar->update(ievent);

            if (! (FR->NTracks()==1))// DA: Check that there is one track per event, or skip
                continue;

            PLTTrack * Track = FR->Track(0);

            //if (Track->Chi2()>12)
            //  continue;

            for (int iroc=0; iroc!= GetNumberOfROCS(telescopeID); iroc++){

                float d_LX = Track->LResidualX(iroc);// DA: local distance between track predicted position and cluster position.
                float d_LY = Track->LResidualY(iroc);

                float cl_LX = -999;
                float cl_LY = -999;

                for (uint16_t icl=0; icl != Track->NClusters(); icl++){
                    if (Track->Cluster(icl)->ROC() == iroc){

                        cl_LX = Track->Cluster(icl)->LX();// DA: local cluster position in the roc
                        cl_LY = Track->Cluster(icl)->LY();


                    }
                }

                // Hits instead of Clusters for Alignment
                // float h_LX = -999;
                // float h_LY = -999;
                // float max_charge = 0;
                // for (int i=0; i != FR->Plane(iroc)->Cluster(0)->NHits(); ++i){
                //
                //   PLTHit * Hit = FR->Plane(iroc)->Cluster(0)->Hit(i);
                //
                //   if (Hit->Charge() > max_charge){
                //     max_charge = Hit->Charge();
                //     h_LX = Hit->LX();
                //     h_LY = Hit->LY();
                //   }
                // }
                //
                // if (fabs(h_LX)>10 || fabs(h_LY)>10)
                //     continue
                //
                // float track_TX = Track->TX(iroc);
                // float track_TY = Track->TY(iroc);
                //
                // float track_LX = Alignment.TtoLX( track_TX, track_TY, 1, iroc);
                // float track_LY = Alignment.TtoLY( track_TX, track_TY, 1, iroc);
                //
                // float d_LX =  (track_LX - h_LX);
                // float d_LY =  (track_LY - h_LY);


                // dX vs dY
                hResidual[iroc].Fill( d_LX, d_LY);
                // DA: take into account only events whose local residuals are less than 1000mm in x and y
                if ((fabs(d_LX) < 10000) && (fabs(d_LY) < 10000)){// DA: Before was < 1000
                    // X vs dY
                    hResidualXdY[iroc].Fill( cl_LX, d_LY);// DA: fill with cluster local position vs residual

                    // Y vs dX
                    hResidualYdX[iroc].Fill( cl_LY, d_LX);

                    gResidualXdY[iroc].SetPoint(gResidualXdY[iroc].GetN(), cl_LX, d_LY );
                    gResidualYdX[iroc].SetPoint(gResidualYdX[iroc].GetN(), cl_LY, d_LX);
                }



            }


        } // end event loop

        float total_angle = 0;
        float total_res = 0;

        for (int iroc=1; iroc!=GetNumberOfROCS(telescopeID); iroc++){
            std::cout << "\nRESIDUALS: " << hResidual[iroc].GetMean(1) << " " << hResidual[iroc].GetMean(2) << std::endl;
            std::cout << "RESIDUALS RMS: " << hResidual[iroc].GetRMS(1) << " " << hResidual[iroc].GetRMS(2) <<std::endl;

            FR->GetAlignment()->AddToLX(1, iroc, hResidual[iroc].GetMean(1));
            FR->GetAlignment()->AddToLY(1, iroc, hResidual[iroc].GetMean(2));

            float angle = atan(hResidualXdY[iroc].GetCorrelationFactor()) ;


            TF1 linear_fun = TF1("","[0]+[1]*x");
            TF1 linear_fun2 = TF1(" ","[0]+[1]*x");
//            gResidualXdY[iroc].Fit(&linear_fun);
//            gResidualYdX[iroc].Fit(&linear_fun2);
            hResidualXdY[iroc].Fit(&linear_fun, "Q");
            hResidualYdX[iroc].Fit(&linear_fun2, "Q");

            float other_angle = atan(linear_fun.GetParameter(1));
            float other_angle2 = atan(linear_fun2.GetParameter(1));
            total_angle += fabs(other_angle);
            total_res += fabs(hResidualXdY[iroc].GetMean(2)) + fabs(hResidualYdX[iroc].GetMean(2));
//            std::cout << "BLA: ANGLE BEFORE: " << FR->GetAlignment()->GetCP(1,iroc)->LR << std::endl;
            FR->GetAlignment()->AddToLR(1, iroc, other_angle);// DA: this was ... other_angle/3.
//            std::cout << "BLA: ANGLE AFTER: " << FR->GetAlignment()->GetCP(1,iroc)->LR << std::endl;

            std::cout << "ROC: " << iroc << " Angle: " << angle << " Other Angle:" << other_angle << " Other Angle 2:" << other_angle2 << std::endl;

//            for (uint8_t i = 0; i != GetNumberOfROCS(telescopeID); i++){
//                printf("%2i   %1i        %15E                       %15E  %15E  %15E\n", 1, i,
//                       FR->GetAlignment()->LR(1,i),
//                       FR->GetAlignment()->LX(1,i),
//                       FR->GetAlignment()->LY(1,i),
//                       FR->GetAlignment()->LZ(1,i) );
//            }


            TCanvas Can;
            Can.cd();

            gResidualXdY[iroc].Draw("AP*");
            Can.SaveAs( OutDir+"/"+TString::Format("gRes%i",iroc) + ".png");

            // 2D Residuals
            hResidual[iroc].SetContour(1024);
            hResidual[iroc].Draw("colz");
            Can.SaveAs( OutDir+"/"+TString(hResidual[iroc].GetName()) + ".png");

            // Residual X-Projection
            gStyle->SetOptStat(1111);
            hResidual[iroc].ProjectionX()->Draw();
            Can.SaveAs( OutDir+"/"+TString(hResidual[iroc].GetName()) + "_X.png");

            // Residual Y-Projection
            hResidual[iroc].ProjectionY()->Draw();
            Can.SaveAs( OutDir+"/"+TString(hResidual[iroc].GetName()) + "_Y.png");

            // 2D Residuals X/dY
            hResidualXdY[iroc].SetContour(1024);
            hResidualXdY[iroc].Draw("colz");
            Can.SaveAs( OutDir+"/"+TString(hResidualXdY[iroc].GetName()) + ".png");

            // 2D Residuals Y/dX
            hResidualYdX[iroc].SetContour(1024);
            hResidualYdX[iroc].Draw("colz");
            Can.SaveAs( OutDir+"/"+TString(hResidualYdX[iroc].GetName()) + ".png");

        }

        std::cout << "END ITERATION " << ialign << " OUT OF 14" << std::endl;

        std::cout << "Sum of magnitude of angle correction per roc: " << total_angle << std::endl;
        std::cout << "Sum of magnitude of residuals in YdX and XdY per rod: " << total_res << std::endl;

        if (total_angle < angle_threshold && total_res < tot_res_threshold){
            std::cout << "total_angle is below " << angle_threshold << " and total_res is below "<< tot_res_threshold << "=> stopping alignment." << std::endl;
            break;
        }



    } // end alignment loop

    std::string outFileName = "NewAlignment.dat";
    std::cout << "saving alignment file \"" << outFileName << "\" with the following parameters:" << std::endl;

    for (uint8_t i = 0; i != GetNumberOfROCS(telescopeID); i++){
        printf("%2i   %1i        %15E                       %15E  %15E  %15E\n", 1, i,
               FR->GetAlignment()->LR(1,i),
               FR->GetAlignment()->LX(1,i),
               FR->GetAlignment()->LY(1,i),
               FR->GetAlignment()->LZ(1,i) );
    }

    FR->GetAlignment()->WriteAlignmentFile(outFileName, FR->NMAXROCS);

    delete FR;

    gROOT->ProcessLine("gErrorIgnoreLevel = 0;");
    return 0;
}
