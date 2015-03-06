#!/usr/bin/env python

"""
Analyze multiple runs from PSI testbeam in May 2014.
Extract the cluster size, leading and second highest charge of  hit as a function of flux.
"""

# ##############################
# Imports
# ##############################

import sys
import math

import ROOT

import RunInfos
from root_style import root_style


###############################
# Read telescope from the commandline
###############################

if not len(sys.argv)==2:
    print "Wrong number of input arguments!"
    print "Usage: {0} telescopeID".format(sys.argv[0])
    sys.exit()
else:
    telescope = int(sys.argv[1])

###############################
# Configuration
###############################

try:
    di_runs = RunInfos.di_di_runs[telescope]
    li_runs_up = RunInfos.di_li_runs_up[telescope]
    li_runs_down = RunInfos.di_li_runs_down[telescope]
    li_runs_down_extended = RunInfos.di_li_runs_down_extended[telescope]
    li_runs_final = RunInfos.di_li_runs_final[telescope]
    di_rocs = RunInfos.di_di_rocs[telescope]
except KeyError:
    print "Invalid telescope! Exiting.."
    sys.exit()


###############################
# Helper: FWHM
###############################

def get_fwhm(h_input):

    h = h_input.Clone()

    debug_this = False

    # First find the maximum bin
    bin_max = h.GetMaximumBin()    
    maximum = h.GetBinContent(bin_max)

    if debug_this:
        print "Maximum {0} found at {1} (=bin: {2})".format(maximum,
                                                            h.GetBinCenter(bin_max),
                                                            bin_max)

    # Go leftwards from the maximum
    # look for the first bin with a value < maximum and record its right edge (technically, the up-edge of the next bin)
    # if none is found, stop at the left edge of the plot    
    for bin_i in reversed(range(1, bin_max)):
        if h.GetBinContent(bin_i) < maximum / 2.:
            break
    bin_left = bin_i    

    left_side_value = h.GetBinLowEdge(bin_left+1)

    if debug_this:
        print "Left edge at at {0} (=bin: {1})".format(left_side_value, bin_left)

    # Go rightwards from the maximum
    # look for the first bin with a value < maximum and record its left edge
    # if none is found, stop at the right edge of the plot    
    for bin_i in range(bin_max+1, h.GetNbinsX()+1):
        if h.GetBinContent(bin_i) < maximum / 2.:
            break
    bin_right = bin_i    

    right_side_value = h.GetBinLowEdge(bin_right)        

    if debug_this:
        print "Right edge at at {0} (=bin: {1})".format(right_side_value, bin_right)

    fwhm = right_side_value - left_side_value

    if debug_this:
        print "FWHM = {0}".format(fwhm)
    
    return fwhm



###############################
# Extract histograms from files
###############################

# Dictionary of lists
# dict_keys = run numbers
# list entries: one per ROC
charge = {}
sum_charge = {}
cluster_size = {}
cluster_size_fraction_1 = {}
cluster_size_fraction_12 = {}
second_charge = {}

#print "{0:>4} {1:>4} {2:>7} {3:>6} {4:>6}".format("Run", "ROC", "Mean", "FWHM", "Entr.", "Gaus", "Left", "Right")
print "{0:>4}   {1:>10}   {2:>4}   {3:>7} {4:>6} {5:>6} {6:>6}".format("Run", "Rate (kHz)", "I", "Mean", "RMS", "Entr.", "SE")


means = {
    "first_three" : [],
    "low_rate"    : [],
    "high_rate"   : [],
}


# Loop over runs
for i_run, run in enumerate(di_runs):

    charge[run] = []
    sum_charge[run] = []
    cluster_size[run] = []
    cluster_size_fraction_1[run] = []
    cluster_size_fraction_12[run] = []
    second_charge[run] = []

    input_rootfile_name = "../plots/000" + str(run) + "/histos.root"
    f = ROOT.TFile.Open(input_rootfile_name)

    for i_roc in range(1, 5):
        # Get the mean cluster size
        h_cs = f.Get("ClusterSize_ROC" + str(i_roc))
        h_cs_proj = h_cs.Project3D("z")
        cluster_size[run].append(h_cs_proj.GetMean())

        # Get the the fraction of events where the cluster size:
        #   is equal to 1
        #   is equal to 1 or 2
        bin_1 = h_cs_proj.FindBin(1)
        bin_2 = h_cs_proj.FindBin(2)
        cluster_size_fraction_1[run].append(1. * h_cs_proj.GetBinContent(bin_1) / h_cs_proj.Integral())
        cluster_size_fraction_12[run].append(1. * h_cs_proj.Integral(bin_1, bin_2) / h_cs_proj.Integral())

        # Get the charge of the second in the cluster
        h_se = f.Get("2ndCharge4_ROC" + str(i_roc))
        second_charge[run].append(h_se.Project3D("z").GetMean())

        # Use the leading charge within a 2-pixel radius
        h_charge = f.Get("1stCharge4_ROC" + str(i_roc))
        charge[run].append(h_charge.Project3D("z").GetMean())

        # Use the sum of  charges within a 4-pixel radius
        h_sumcharge = f.Get("SumCharge4_ROC" + str(i_roc))

        htmp = h_sumcharge.Project3D("z")
        htmp.SetBinContent(1, 0)
        sum_charge[run].append(htmp.GetMean())
        

        
        
        rs = root_style()
        rs.set_style(1000,1000,1)
        ROOT.gROOT.ForceStyle()
        
        c = rs.get_canvas("")
            
        #if run in [322, 325, 327, 330, 333]:
        if True:
            h_sumcharge_proj_z = h_sumcharge.Project3D("z")
            if i_roc in [4]:
                
                for i_iter in range(-1,3):

                    # Ignore the zero-bin
                    if i_iter == 0:
                        h_sumcharge_proj_z.SetBinContent(1,0)                    


                        if run in [322, 325, 327]:
                            means["first_three"].append(h_sumcharge_proj_z.GetMean()) 
                        if run in [322, 347, 348]:
                            means["low_rate"].append(h_sumcharge_proj_z.GetMean())
                        if run in [333, 350, 352]:
                            means["high_rate"].append(h_sumcharge_proj_z.GetMean())
                           

                    if i_iter == 1:
                        # Highest Bin
                        h_sumcharge_proj_z.SetBinContent( h_sumcharge_proj_z.GetNbinsX(), 0)
                        # Overflow Bin
                        h_sumcharge_proj_z.SetBinContent( h_sumcharge_proj_z.GetNbinsX()+1, 0)

                    if i_iter == 2:
                        # Bins below 7500
                        for i_bin in range( h_sumcharge_proj_z.FindBin(7500)+1):
                            h_sumcharge_proj_z.SetBinContent(i_bin, 0)
                    

                    entries = h_sumcharge_proj_z.Integral()
                    mean = h_sumcharge_proj_z.GetMean()
                    rms  = h_sumcharge_proj_z.GetRMS()
                    #fwhm = get_fwhm(h_sumcharge_proj_z)
                    #sigma = fwhm/2.
                    #
                    #bin_4_sigma_left = h_sumcharge_proj_z.FindBin(mean - 4*sigma)
                    #bin_3_sigma_left = h_sumcharge_proj_z.FindBin(mean - 3*sigma)
                    #bin_3_sigma_right = h_sumcharge_proj_z.FindBin(mean + 3*sigma)
                    #bin_4_sigma_right = h_sumcharge_proj_z.FindBin(mean + 4*sigma)
                    #
                    #integral_left = h_sumcharge_proj_z.Integral(bin_4_sigma_left, bin_3_sigma_left)
                    #integral_right = h_sumcharge_proj_z.Integral(bin_3_sigma_right, bin_4_sigma_right)
                    #
                    #f_gauss_3_to_4 = (0.2699796-0.006334)/(2*100)
                    #print "{0:>4} {1:>4} {2:>7.1f} {3:>6.0f} {4:>6.0f} {5:>6.1f} {6:>6.1f} {7:>6.1f}".format(run, i_roc, mean, fwhm, entries, entries * f_gauss_3_to_4, integral_left, integral_right)

                    print "{0:>4}   {1:>10.0f}   {2:>4}   {3:>7.1f} {4:>6.0f} {5:>6.0f} {6:>6.0f}".format(run, di_runs[run], i_iter, mean, rms, entries, rms/math.sqrt(entries))

                

    # End of loop over ROCs
# End loop over runs

print means
print "{0:>15}: {1:>10} {2:>10}".format("Runs", "Mean", "RMS")
for k,v in means.iteritems():

    m = 1. * sum(v)/len(v)
    rms = math.sqrt(1.*sum([pow(x-m,2) for x in v])/len(v))
    
    print "{0:>15}: {1:>10.1f} {2:>10.1f}".format(k, m, rms)


###############################
# Exceptions: Variable Error
###############################

class VariableError(Exception):
    def __str__(self):
        return "Do not know how to handle this plotting variable."


###############################
# make_plots
###############################

def make_plots():

    # Prepare pretty ROOT
    rs = root_style()
    rs.set_style(1000,1000,1)
    ROOT.gROOT.ForceStyle()

    c = rs.get_canvas("")
    c.SetLogx(1)



    for plot_var in ["cluster_size", 
                     #"cluster_size_fraction_1", 
                     #"cluster_size_fraction_12", 
                     "charge", 
                     "sum_charge", 
                     #"second_charge"
    ]:

        if plot_var == "cluster_size":
            di_values = cluster_size
            legend_origin_x = 0.3
            legend_origin_y = 0.4
        elif plot_var == "cluster_size_fraction_1":
            di_values = cluster_size_fraction_1
            legend_origin_x = 0.3
            legend_origin_y = 0.2
        elif plot_var == "cluster_size_fraction_12":
            di_values = cluster_size_fraction_12
            legend_origin_x = 0.3
            legend_origin_y = 0.2
        elif plot_var == "charge":
            di_values = charge
            legend_origin_x = 0.3
            legend_origin_y = 0.8
        elif plot_var == "sum_charge":
            di_values = sum_charge
            legend_origin_x = 0.3
            legend_origin_y = 0.85
        elif plot_var == "second_charge":
            di_values = second_charge
            legend_origin_x = 0.3
            legend_origin_y = 0.65
        else:
            raise VariableError(plot_var)

        legend_size_x = 0.1
        legend_size_y = 0.045 * 3




        # Loop over ROCs
        for i_roc in range(1, 5):

            if plot_var == "sum_charge" and i_roc ==4:
                legend_origin_y = 0.5
            
            # Prepare Legend
            legend = rs.make_legend(legend_origin_x, 
                                    legend_origin_y,
                                    3)
                                
            if plot_var == "cluster_size":
                h = ROOT.TH2F("", "", 100, 1, 40000, 100, 0., 3)
                h.GetYaxis().SetTitle("Mean cluster size")
            elif plot_var == "cluster_size_fraction_1":
                h = ROOT.TH2F("", "", 100, 1, 40000, 100, 0., 1.)
                h.GetYaxis().SetTitle("Fraction of size == 1 clusters")
            elif plot_var == "cluster_size_fraction_12":
                h = ROOT.TH2F("", "", 100, 1, 40000, 100, 0., 1.)
                h.GetYaxis().SetTitle("Fraction of 1 <= size <= 2 clusters")
            elif plot_var == "charge":
                h = ROOT.TH2F("", "", 100, 1, 40000, 100, 0., 30000)
                h.GetYaxis().SetTitle("Leading charge within 4-pixel radius")
            elif plot_var == "sum_charge":
                h = ROOT.TH2F("", "", 100, 1, 40000, 100, 0., 30000)
                h.GetYaxis().SetTitle("Sum of charge within 4-pixel radius")
            elif plot_var == "second_charge":
                h = ROOT.TH2F("", "", 100, 1, 40000, 100, 0., 15000)
                h.GetYaxis().SetTitle("Charge of second hit in cluster")
            else:
                raise VariableError

            h.GetXaxis().SetTitle("Flux [kHz/cm^{2}]")
            #

            #h.GetXaxis().SetTitleSize(0.06)
            #h.GetYaxis().SetTitleSize(0.06)
            #h.GetXaxis().SetLabelSize(0.06)
            #h.GetYaxis().SetLabelSize(0.06)

            h.Draw()

            # Prepare efficiency TGraphs
            li_grs = []

            if plot_var ==  "cluster_size":
                pass
            else:
                h.GetYaxis().SetTitleOffset(1.95)



            for direction in ["up"]:

                # Choose runs to use
                if direction == "up":
                    li_runs = li_runs_up
                elif direction == "down":
                    li_runs = li_runs_down
                else:
                    li_runs = li_runs_final

                # Initialize TGraph objects. For tracking we need more points for multiple event slices
                gr = ROOT.TGraphErrors(len(li_runs))

                # Loop over runs
                for irun, run in enumerate(sorted(li_runs)):
                    gr.SetPoint(irun, di_runs[run], di_values[run][i_roc-1])  # i_roc-1 since we have no data for ROC0
                    gr.SetPointError(irun, 0, 490) # two sigma

                # Make things look nice:
                # Legend Entries

                if direction == "up":
                    legend.AddEntry(gr, "Increasing flux", "P")
                elif direction == "down":
                    legend.AddEntry(gr, "Decreasing flux", "P")
                elif direction == "final":
                    legend.AddEntry(gr, "Highest flux", "P")

                # Markers
                # going up
                if direction == "up":
                    gr.SetMarkerStyle(22)
                    gr.SetMarkerColor(ROOT.kRed)
                #  down
                elif direction == "down":
                    gr.SetMarkerStyle(23)
                    gr.SetMarkerColor(ROOT.kBlue)
                # final high flux
                else:
                    gr.SetMarkerStyle(21)
                    gr.SetMarkerColor(ROOT.kGreen)

                # Protect graphs from autodelete
                li_grs.append(gr)
                gr.SetMarkerSize(3)
                gr.Draw("PSAME")
                legend.Draw()

                if plot_var == "cluster_size":
                    outfile_name = "ClusterSize"
                elif plot_var == "cluster_size_fraction_1":
                    outfile_name = "ClusterSizeFraction1"
                elif plot_var == "cluster_size_fraction_12":
                    outfile_name = "ClusterSizeFraction12"
                elif plot_var == "charge":
                    outfile_name = "LeadingCharge"
                elif plot_var == "sum_charge":
                    outfile_name = "SumCharge"
                elif plot_var == "second_charge":
                    outfile_name = "SecondCharge"
                else:
                    raise VariableError
        
            outfile_name += "_Telescope{0}_{1}_{2}".format(telescope, di_rocs[i_roc], "all")

            c.Print(outfile_name + ".png")
            c.Print(outfile_name + ".pdf")

# End of Make Plots

make_plots()
