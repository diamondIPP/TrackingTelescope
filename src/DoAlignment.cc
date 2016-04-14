#include "DoAlignment.h"

using namespace std;

/** DoAlignment: Produce alignment constants and save them to NewAlignment.dat */
int DoAlignment (std::string const InFileName,
                 TFile * out_f,
                 TString const RunNumber,
                 int telescopeID)
{
    TString const PlotsDir = "plots/";
    TString const OutDir = PlotsDir + RunNumber;

    gStyle->SetOptStat(0);

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
                                   GetNumberOfROCS(telescopeID), GetUseGainInterpolator(telescopeID), GetUseExternalCalibrationFunction(telescopeID), true);
    }
    else{
        FR = new PSIBinaryFileReader(InFileName, GetCalibrationFilename(telescopeID), GetAlignmentFilename(telescopeID, true),
                                     GetNumberOfROCS(telescopeID), GetUseGainInterpolator(telescopeID), GetUseExternalCalibrationFunction(telescopeID) );
        ((PSIBinaryFileReader*) FR)->CalculateLevels(10000, OutDir);
    }

    uint32_t stopAt = ((PSIRootFileReader*) FR)->fTree->GetEntries();
//    uint32_t stopAt = 5e4;

    /** Apply Masking */
    FR->ReadPixelMask(GetMaskingFilename(telescopeID));

    for (int ialign=0; ialign!=2; ialign++){

        for (uint8_t iroc=1; iroc!=GetNumberOfROCS(telescopeID) - 1; iroc++){
            FR->GetAlignment()->AddToLX( 1, iroc, x_align[iroc] );
            FR->GetAlignment()->AddToLY( 1, iroc, y_align[iroc] );
            //FR->GetAlignment()->AddToLR( 1, iroc, r_align[iroc] );
        }

        for (int iroc_align = 1; iroc_align != GetNumberOfROCS(telescopeID); ++iroc_align) {

            std::cout << "GOING TO ALIGN: " << iroc_align << std::endl;

            FR->ResetFile();
            FR->SetPlaneUnderTest( iroc_align );

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
                hResidual.push_back(    TH2F(  Form("Residual_ROC%i",iroc),     Form("Residual_ROC%i",iroc), 200, -.2, .2, 200, -.2, .2));
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

                /** loop over hits */
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

                float track_LX = FR->GetAlignment()->TtoLX( track_TX, track_TY, 1, iroc_align);
                float track_LY = FR->GetAlignment()->TtoLY( track_TX, track_TY, 1, iroc_align);

                float d_LX =  (track_LX - h_LX);
                float d_LY =  (track_LY - h_LY);

                //std::cout << "Track LX/LY" << track_LX << " " << track_LY << std::endl;

                if ( !(fabs(d_LX)<2) || !(fabs(d_LY)<2) ) continue;

                hResidual[iroc_align].Fill( d_LX, d_LY); // dX vs dY

                hResidualXdY[iroc_align].Fill( h_LX, d_LY); // X vs dY

                hResidualYdX[iroc_align].Fill( h_LY, d_LX); // Y vs dX

            } // end event loop
            cout << "Loop duration:" << (clock() - now) / CLOCKS_PER_SEC << endl;

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
            hResidual[iroc_align].Draw("colz");
            Can.SaveAs( OutDir+"/"+TString(hResidual[iroc_align].GetName()) + ".gif");

            // Residual X-Projection
            gStyle->SetOptStat(1111);
            hResidual[iroc_align].ProjectionX()->Draw();
            Can.SaveAs( OutDir+"/"+TString(hResidual[iroc_align].GetName()) + "_X.gif");

            // Residual Y-Projection
            hResidual[iroc_align].ProjectionY()->Draw();
            Can.SaveAs( OutDir+"/"+TString(hResidual[iroc_align].GetName()) + "_Y.gif");

            // 2D Residuals X/dY
            hResidualXdY[iroc_align].Draw("colz");
            Can.SaveAs( OutDir+"/"+TString(hResidualXdY[iroc_align].GetName()) + ".gif");

            // 2D Residuals Y/dX
            hResidualYdX[iroc_align].Draw("colz");
            Can.SaveAs( OutDir+"/"+TString(hResidualYdX[iroc_align].GetName()) + ".gif");

            for (uint8_t i = 0; i != GetNumberOfROCS(telescopeID); i++){
                std::cout << int(i) << setprecision(7) << " " << x_align[i] << " " << y_align[i] << " " << z_align[i] << " " << r_align[i] <<std::endl;
            }
        } // end loop over rocs

        FR->GetAlignment()->WriteAlignmentFile("NewAlignment.dat", FR->NMAXROCS);

    } // end alignment loop




    std::cout << "PART TWO!!!!!" << std::endl;


    for (int ialign=1; ialign!=15;ialign++){

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
        }

        // Event Loop
        for (uint32_t ievent = 0; FR->GetNextEvent() >= 0; ++ievent) {

            if (ievent > stopAt)
                break;

            /** print progress */
            print_progress(ievent, stopAt);

            if (! (FR->NTracks()==1))
                continue;

            PLTTrack * Track = FR->Track(0);

            //if (Track->Chi2()>12)
            //  continue;

            for (int iroc=0; iroc!= GetNumberOfROCS(telescopeID); iroc++){

                float d_LX = Track->LResidualX(iroc);
                float d_LY = Track->LResidualY(iroc);

                float cl_LX = -999;
                float cl_LY = -999;

                for (uint16_t icl=0; icl != Track->NClusters(); icl++){
                    if (Track->Cluster(icl)->ROC() == iroc){

                        cl_LX = Track->Cluster(icl)->LX();
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

                if ((fabs(d_LX) < 1000) && (fabs(d_LY) < 1000)){
                    // X vs dY
                    hResidualXdY[iroc].Fill( cl_LX, d_LY);

                    // Y vs dX
                    hResidualYdX[iroc].Fill( cl_LY, d_LX);

                    gResidualXdY[iroc].SetPoint(gResidualXdY[iroc].GetN(), cl_LX, d_LY );
                }



            }


        } // end event loop

        for (int iroc=1; iroc!=GetNumberOfROCS(telescopeID); iroc++){
            std::cout << "RESIDUALS: " << hResidual[iroc].GetMean(1) << " " << hResidual[iroc].GetMean(2) << std::endl;
            std::cout << "RESIDUALS RMS: " << hResidual[iroc].GetRMS(1) << " " << hResidual[iroc].GetRMS(2) <<std::endl;

            FR->GetAlignment()->AddToLX(1, iroc, hResidual[iroc].GetMean(1));
            FR->GetAlignment()->AddToLY(1, iroc, hResidual[iroc].GetMean(2));

            float angle = atan(hResidualXdY[iroc].GetCorrelationFactor()) ;


            TF1 linear_fun = TF1("","[0]+[1]*x");
            gResidualXdY[iroc].Fit(&linear_fun);


            float other_angle = atan(linear_fun.GetParameter(1));

            FR->GetAlignment()->AddToLR(1, iroc, other_angle/3.);

            std::cout << "ROC: " << iroc << " Angle: " << angle << " Other Angle:" << other_angle << std::endl;

            for (uint8_t i = 0; i != GetNumberOfROCS(telescopeID); i++){
                printf("%2i   %1i        %15E                       %15E  %15E  %15E\n", 1, i,
                       FR->GetAlignment()->LR(1,i),
                       FR->GetAlignment()->LX(1,i),
                       FR->GetAlignment()->LY(1,i),
                       FR->GetAlignment()->LZ(1,i) );
            }


            TCanvas Can;
            Can.cd();

            gResidualXdY[iroc].Draw("AP*");
            Can.SaveAs( OutDir+"/"+TString::Format("gRes%i",iroc) + ".gif");

            // 2D Residuals
            hResidual[iroc].Draw("colz");
            Can.SaveAs( OutDir+"/"+TString(hResidual[iroc].GetName()) + ".gif");

            // Residual X-Projection
            gStyle->SetOptStat(1111);
            hResidual[iroc].ProjectionX()->Draw();
            Can.SaveAs( OutDir+"/"+TString(hResidual[iroc].GetName()) + "_X.gif");

            // Residual Y-Projection
            hResidual[iroc].ProjectionY()->Draw();
            Can.SaveAs( OutDir+"/"+TString(hResidual[iroc].GetName()) + "_Y.gif");

            // 2D Residuals X/dY
            hResidualXdY[iroc].Draw("colz");
            Can.SaveAs( OutDir+"/"+TString(hResidualXdY[iroc].GetName()) + ".gif");

            // 2D Residuals Y/dX
            hResidualYdX[iroc].Draw("colz");
            Can.SaveAs( OutDir+"/"+TString(hResidualYdX[iroc].GetName()) + ".gif");

        }



    } // end alignment loop

    FR->GetAlignment()->WriteAlignmentFile("NewAlignment.dat", FR->NMAXROCS);

    delete FR;

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
