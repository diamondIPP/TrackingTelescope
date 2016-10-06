#include "DoAlignment.h"

using namespace std;

/** DoAlignment: Produce alignment constants and save them to NewAlignment.dat */
int DoAlignment (std::string const InFileName, TString const RunNumber, int telescopeID)
{
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
            FR->SetPlaneUnderTest( iroc_align );// ignore plane for tracking

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

                if (ievent > stopAt)
                    break;

                /** print progress */
                print_progress(ievent, stopAt);

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

                float track_TX = Track->TX(iroc_align);
                float track_TY = Track->TY(iroc_align);

                float track_LX = FR->GetAlignment()->TtoLX( track_TX, track_TY, 1, iroc_align);// Local position of the track in the plane under test
                float track_LY = FR->GetAlignment()->TtoLY( track_TX, track_TY, 1, iroc_align);

                float d_LX =  (track_LX - h_LX);// residuals of track local position and highes charge hit local position
                float d_LY =  (track_LY - h_LY);

                //std::cout << "Track LX/LY" << track_LX << " " << track_LY << std::endl;
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
            print_progress(ievent, stopAt);

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

void print_progress(uint32_t ievent, uint32_t stopAt){

    if (ievent % 10 == 0 && ievent > 1000){
        if (ievent != 0) cout << "\r";
        cout << "Processed events: "  << setprecision(2) << setw(5) << setfill('0') << fixed << float(ievent) / stopAt * 100 << "% ";
        if ( stopAt - ievent < 10 ) cout << "|" <<string( 50 , '=') << ">";
        else cout << "|" <<string(int(float(ievent) / stopAt * 100) / 2, '=') << ">";
        if ( stopAt - ievent < 10 ) cout << "| 100%    ";
        else cout << string(50 - int(float(ievent) / stopAt * 100) / 2, ' ') << "| 100%    ";
        cout.flush();
        if ( stopAt - ievent < 10 ) cout << endl;
    }
}
