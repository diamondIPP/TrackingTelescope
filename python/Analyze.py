#!/usr/bin/env python

"""
Analyze multiple runs from PSI testbeam in May 2014
"""

###############################
# Imports
###############################

from math import log10, floor
import sys

import ROOT

import RunInfos
from root_style import root_style



###############################
# Helpers
###############################

def round_to_2(x):
    return round(x, -int(floor(log10(x)))+1)


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
    li_runs_up = RunInfos.di_li_runs_up_less[telescope]
    li_runs_down = RunInfos.di_li_runs_down_less[telescope]
    di_rocs = RunInfos.di_di_rocs[telescope]
except KeyError:
    print "Invalid telescope! Exiting.."
    sys.exit()


li_names = []
#li_names +=  ["PulseHeightTrack6_ROC"+str(iroc)+"_All" for iroc in range(6)]
#li_names +=  ["PulseHeightTrack6_ROC"+str(iroc)+"_NPix1" for iroc in range(6)]
#li_names +=  ["PulseHeightTrack6_ROC"+str(iroc)+"_NPix2" for iroc in range(6)]
#li_names +=  ["PulseHeightTrack6_ROC"+str(iroc)+"_NPix3Plus" for iroc in range(6)]


li_names +=  ["1stCharge4_ROC{0}_z".format(iroc) for iroc in range(1,5)]
li_names +=  ["SumCharge4_ROC{0}_z".format(iroc) for iroc in range(1,5)]

###############################
# Extract histograms from files
###############################

di_histos = {}

for name in li_names:
    di_histos[name] = {}

for i_run, run in enumerate(li_runs_up + li_runs_down):

    print "Doing", run
    input_rootfile_name = "../plots/000"+str(run)+"/histos.root"
    f = ROOT.TFile.Open( input_rootfile_name )

    for name in li_names:
        h = f.Get( name ).Clone()
        h.SetDirectory(0)
        if h.Integral()>0:
            h.Scale(100./h.Integral())
        h.SetTitle("")
        h.GetXaxis().SetNdivisions(505)
        if "Charge" in name:
            h.GetXaxis().SetTitle("signal a.u.")
        h.GetYaxis().SetTitle("Fraction of clusters [%]")
        di_histos[name][run]=h


###############################
# Prepare pretty ROOT
###############################

rs = root_style()
rs.set_style(1000,1000,1)
ROOT.gROOT.ForceStyle()

#ROOT.gStyle.SetPadLeftMargin(0.15)
#ROOT.gStyle.SetPadBottomMargin(0.15)
#ROOT.gStyle.SetPadRightMargin(0.05)
#ROOT.gStyle.SetPadTopMargin(0.05)
#ROOT.gROOT.ForceStyle()


c = rs.get_canvas("")



#c.SetGrid(1,1)

li_colors = [ROOT.kRed,      ROOT.kBlue+1,     ROOT.kBlack,
             ROOT.kOrange-1, ROOT.kViolet+1,   ROOT.kGreen+1,
             ROOT.kGray,     ROOT.kYellow,
         ]*10


def GetMaximumExceptBin(h, ibin=1):

  if (h.GetMaximumBin() == ibin):
    return h.GetMaximum(h.GetMaximum())
  else:
    return h.GetMaximum()


###############################
# Combine histos into plot
###############################

for name in li_names:
    for direction in ["up"]:

        # Choose runs to use
        if direction == "up":
            li_runs = li_runs_up
        elif direction == "down":
            li_runs = li_runs_down

        legend_origin_x     = 0.55
        legend_origin_y     = 0.8

        legend = rs.make_legend(legend_origin_x, 
                                legend_origin_y,
                                len(li_runs))
        legend.SetTextAlign(32)
        legend.SetMargin(.4)

        the_max = max( [GetMaximumExceptBin(di_histos[name][run])for run in li_runs])

        print h.GetNbinsX()
        
        x = ROOT.Long()
        h = di_histos[name][li_runs[0]]
        h.GetBinWithContent(GetMaximumExceptBin(h),x)
        f = h.GetXaxis().GetXmax()/h.GetBinCenter(x)
        print f

        #print max_bin

        the_max *= 1.1

        for i_run, run in enumerate(li_runs):
        

            di_histos[name][run].GetXaxis().SetLimits(0, f)
            
        

            di_histos[name][run].SetLineColor( li_colors[i_run])
            rate = int(round((round_to_2(di_runs[run]))))
            legend.AddEntry( di_histos[name][run], "{0} kHz/cm".format(rate)+"^{2}", "L" )

            di_histos[name][run].SetMinimum(0)
            di_histos[name][run].SetMaximum(the_max)
            #di_histos[name][run].GetXaxis().SetRange(0, 49900)


            di_histos[name][run].GetXaxis().SetRangeUser(0, 3)
            
            #di_histos[name][run].GetXaxis().SetTitleSize(0.06)
            #di_histos[name][run].GetYaxis().SetTitleSize(0.06)
            #di_histos[name][run].GetXaxis().SetLabelSize(0.06)
            #di_histos[name][run].GetYaxis().SetLabelSize(0.06)

            if i_run == 0:
                di_histos[name][run].Draw()
            else:
                di_histos[name][run].Draw("SAME")

        legend.Draw()

        # Extract the ROC number from the name
        roc_pos_in_string = name.find("ROC") + 3 # because ROC has three letters
        roc_index = int(name[roc_pos_in_string])
        

        c.Print("plots/{0}_Telescope{1}_{2}_{3}.pdf".format(name, telescope, di_rocs[roc_index],  direction))
