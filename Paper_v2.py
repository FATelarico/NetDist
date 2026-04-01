#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Title: graph-tool SBM

 Purpose : This script runs the SBM with graph-tool for the
           data on FDI flows

Project script number: #0001

Author: Fabio A. Telarico
Contact details: Fabio-Ashtar.Telarico@fdvsuni-lj.si
Date script created: 2026-01-19 3:47 pm CET
Date script last modified: 2026-##-## ##:## ## CET
"""

# %% Setup

# %%% Libraries
import os
import __main__
from graph_tool.all import *
import re
import matplotlib
import numpy
import pickle
import json
import pathlib


# %%% Query working directory
if(hasattr(__main__, "__file__")):
    wd = os.path.abspath(__file__) # Get script's location while in CLI
    wd = os.path.dirname(wd) # Extract the folder
    is_interactive = False
    # print('In an shell')
    # print(wd) # Check
else:
    is_interactive = True
    print('In an interactive environment: SET THE WD MANUALLY!')
    wd = "/home/ashtar/JIBE_AfricaNewGlobalOrder/" # Dummy value
    
if os.getcwd()!= wd:
    os.chdir(wd)
    print(os.getcwd())

# %%% Parellisation
# os.environ["OMP_WAIT_POLICY"] = "active"
# os.environ["OMP_WAIT_POLICY"] = "passive" # Default!
# #| Note: Due to an OpenMP API limitation, this CANNOT be changed
# #| after graph-tool has been imported.

if(graph_tool.openmp_get_num_threads()==1):
    graph_tool.openmp_set_num_threads(64)
    
if(graph_tool.openmp_get_thresh()>300): # Restore default value
    graph_tool.openmp_set_thresh(300)

# %% 1. Import graph from R
if(is_interactive):
    data = load_graph('./Input/lnBlss.gml',
                      fmt = "graphml", ignore_vp=['_graphml_vertex_id'],
                      ignore_ep = ['_graphml_edge_id'])
    
    data.list_properties() # Check that `name`, `t`, and `weight` are present

# %% 2. Blockmodel without DC

# %%% 2.1 Run BM
if(is_interactive):
    seed_rng(74651) # Set a seed for the RNG
    
    # Estimate hierarchical blockmodel edge weights and degree correction
    res_woutDC = minimize_nested_blockmodel_dl(data, state_args=dict(
        deg_corr = False, recs = [data.ep.weight], rec_types = ['real-normal'],
        pclabel = data.vp.t
    ))
    
    with open('./Output/res_HBSBM-woutDC.pkl', 'wb') as file:
        # Save to file
        pickler = pickle.Pickler(file)
        pickler.dump(res_woutDC)
         

# %%% 2.2 Improve solution with merge-split
if os.path.exists('./Output/res_HSBM-woutDC-MCMC.pkl'):
    print('Nothing to do without DC')
else:
    if os.path.exists('./Output/res_HSBM-woutDC-temp.pkl'):
        with open('./Output/res_woutDC-temp.pkl', "rb") as file:
            res_woutDC = pickle.load(file)
        # i_start =  # Run 1 on 2026/01/xx
        print(f"Temporary result loaded ()")
    else:
        with open('./Output/res_HBSBM-woutDC.pkl', 'rb') as file:
            res_woutDC = pickle.load(file)
        i_start = 0
        print("Base result loaded")
    
    MDL = res_woutDC.entropy()
    print(f"Starting from {i_start}")
    for i in range(i_start, 10000): 
        seed_rng(49865+i) # Set a seed for the RNG
        res_woutDC.multiflip_mcmc_sweep(niter=100)
        if(res_woutDC.entropy()<MDL):
            print(f"Sweep {i}: Description length = {res_woutDC.entropy():.4f}")
            MDL = res_woutDC.entropy()
        if i % 100 == 0:
            with open('./Output/res_woutDC-temp.pkl', 'wb') as file:
                pickler = pickle.Pickler(file)
                pickler.dump(res_woutDC)
        
    # Save to file
    with open('./Output/res_HSBM-woutDC-MCMC.pkl', 'wb') as file:
        pickler = pickle.Pickler(file)
        pickler.dump(res_woutDC)

# %% 3. Blockmodel with DC

# %%% 3.1 Run BM
if(is_interactive):
    seed_rng(74651) # Set a seed for the RNG
    
    # Estimate hierarchical blockmodel edge weights and degree correction
    res_withDC = minimize_nested_blockmodel_dl(data, state_args=dict(
        deg_corr = True, recs = [data.ep.weight], rec_types = ['real-normal'],
        pclabel = data.vp.t
    ))
    
    with open('./Output/res_HBSBM-withDC.pkl', 'wb') as file:
        # Save to file
        pickler = pickle.Pickler(file)
        pickler.dump(res_withDC)
         

# %%% 3.2 Improve solution with merge-split
if os.path.exists('./Output/res_HSBM-withDC-MCMC.pkl'): # If already optimised
    print('Nothing to run!')
else:
    if os.path.exists('./Output/res_withDC-temp.pkl'): # If temporary file
        # Load it
        with open('./Output/res_withDC-temp.pkl', "rb") as file:
            res_withDC = pickle.load(file)
        i_start = pathlib.Path('./Output/res_withDC-idx.txt')
        i_start = i_start.read_text(encoding="utf-8").splitlines()
        i_start = int(next(ii for ii in reversed(i_start) if ii.strip()))
        print(f"Temporary result loaded (i={i_start})")
    else:
        # Withou temporary file, laod the main result
        with open('./Output/res_HBSBM-withDC.pkl', 'rb') as file:
            res_withDC = pickle.load(file)
        i_start = 0
        print("Base result loaded")
    for i in range(i_start, 10000): # Iterations for improvement
        d = res_withDC.entropy()
        seed_rng(132564+i) # Set a seed for the RNG
        res_withDC.multiflip_mcmc_sweep(niter=1000)
        if(res_withDC.entropy()<d):
            print(f"Sweep {i}: New MDL = {res_withDC.entropy():.4f}")
        if i % 100 == 0: # Every one hundred iterations
            # Write down counter
            with open('./Output/res_withDC-idx.txt', 'a', encoding='utf-8') as f:
                f.write(f"{i}\n")
            # Export temporary result
            with open('./Output/res_withDC-temp.pkl', 'wb') as file:
                pickler = pickle.Pickler(file)
                pickler.dump(res_withDC)
    # Save optimised result to file
    with open('./Output/res_HSBM-withDC-MCMC.pkl', 'wb') as file:
        pickler = pickle.Pickler(file)
        pickler.dump(res_withDC)
