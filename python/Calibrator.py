#!/usr/bin/env python

"""
Redo the GainCal fit
"""

###############################
# Imports
###############################

import sys

import ROOT


###############################
# Configuration
###############################

infile_name = sys.argv[1] #"../Calibrations/phCalibration_C1.dat"
outfile_name = infile_name.replace("phCalibration","phCalibrationGErfFit")

lo_range_string = "Low range:"
hi_range_string = "High range:"

hi_range_factor = 7

produce_debugplots = True

output_header = "Parameters of the vcal vs. pulse height fits"

fit_fun_string = "[3]*(TMath::Erf((x-[0])/[1])+[2])"


###############################
# Write the output header
###############################

out_f = open(outfile_name, "w")

out_f.write(output_header + "\n")
out_f.write(fit_fun_string + "\n\n")


###############################
# Read the header
###############################

in_f = open(infile_name, "r")


# Ignore first line
in_f.readline()

# Second line: low range
line = in_f.readline()

# make sure it is the low range line
if not lo_range_string in line:
    print "Error! Line is not a low-range line: ", line
    sys.exit()

# clean the line
line = line.replace(lo_range_string, "").strip()
lo_ranges = [int(x) for x in line.split(" ")]

# Third line: high range
line = in_f.readline()

# make sure it is the high range line
if not hi_range_string in line:
    print "Error! Line is not a high-range line: ", line
    sys.exit()

# clean the line
line = line.replace(hi_range_string, "").strip()
hi_ranges = [int(x) for x in line.split(" ")]

# high range gets an extra factor
hi_ranges = [x*7 for x in hi_ranges]

vcal_values = lo_ranges + hi_ranges

print "Got vcal_values: ", vcal_values

# Ignore the empty line
in_f.readline()

c1 = ROOT.TCanvas("","", 800, 800)

for line in in_f:
    
    # Format of line should be:
    # -128 -104 -78 -53 -24 -123 -95 -72 -50  32 104 150 166 169    Pix  0  0
    # So we get the index of the "Pix" position
    # and make sure we have as many data entries as vcal values
    
    # get the line and split using spaces, remove empty fields
    line = line.strip()    
    atoms = line.split(" ")
    atoms = [x for x in atoms if not x==""]
    
    # make sure we have the right number of data entries
    pix_index = atoms.index("Pix")
    if not  pix_index == len(vcal_values):
        print "Mismatch data/vcal_values in line:", line
        sys.exit()
        
    gr = ROOT.TGraph()
    
    for i in range(pix_index):
        if not atoms[i] == "0":
            gr.SetPoint(gr.GetN(), vcal_values[i], float(atoms[i]))

    fun = ROOT.TF1("fun", fit_fun_string)
    fun.SetParameter(0, 300)
    fun.SetParameter(1, 350)
    fun.SetParameter(2,-0.5)
    fun.SetParameter(3, 150)
    
    gr.Fit(fun)
    
    if produce_debugplots:
        gr.Draw("APL*")    
        c1.Print("debugplots/Pix_{0}_{1}.pdf".format(atoms[pix_index+1], atoms[pix_index+2]))

    out_string = "{0} {1} {2} {3}    Pix {4} {5}".format( fun.GetParameter(0),
                                                          fun.GetParameter(1),
                                                          fun.GetParameter(2),
                                                          fun.GetParameter(3),
                                                          atoms[pix_index+1], 
                                                          atoms[pix_index+2])
    out_f.write(out_string + "\n")
