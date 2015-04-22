#!/usr/bin/env python

"""
Create a ROOT file for testing a ROOT file reader in the analysis
"""

###############################
# Imports
###############################

import time
import array
import random

import ROOT


###############################
# Configuration
###############################

outfile_name = "test.root"
n_events = 1000

scalar_branches = [
    ["event_number", int],
    ["time", float],
]

vector_branches = [
    ["plane",  int],
    ["col",    int],
    ["row",    int],
    ["adc",    int],
    ["charge", int]
]    


###############################
# Prepare Tree
###############################

# Store all branches here
d = {}

# Create file/tree
outfile = ROOT.TFile( outfile_name, 'recreate')
tree = ROOT.TTree("tree", "tree")

# Add scalar branches                   
for branch in scalar_branches:
    name = branch[0]
    datatype = branch[1]
    
    if datatype == int:
        d[name] = array.array('i', [0])
        setattr(tree, name, d[name])        
        tree.Branch(name, d[name], name+"/I")          
    elif datatype == float:
        d[name] = array.array('f', [0.])
        setattr(tree, name, d[name])        
        tree.Branch(name, d[name], name+"/F")          
# End of adding scalar branches

# Add vector branches
for branch in vector_branches:
    name = branch[0]
    datatype = branch[1]

    d[name] = ROOT.std.vector( datatype )()
    setattr(tree, name, d[name])
    tree.Branch(name, d[name])          
# End of adding vector branches    


###############################
# Fill Tree
###############################
    
for i in xrange(n_events):

    # Reset branches
    for b in scalar_branches:
        name = b[0]
        d[name][0] = 0
    for b in vector_branches:
        name = b[0]
        d[name].resize(0)

    # Fill per event branches
    d["event_number"][0] = i
    d["time"][0] = time.clock()

    # Fill per hit branches    
    # first decide how many hits we want
    for i_plane in range(4):
        d["plane"].push_back(i_plane)
        d["col"].push_back( random.randint(20,50))
        d["row"].push_back( random.randint(20,50))
        d["adc"].push_back( random.randint(80,120))
        d["charge"].push_back( random.randint(3000,20000))
    # End of per-hit loop

    tree.Fill()    

# End of Event Loop

# Save everything & exit cleanly
tree.AutoSave()
outfile.Close()    
