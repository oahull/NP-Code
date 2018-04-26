#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 17:22:23 2017

@author: oah
"""

def parse_files(prefix, file_range):
    '''
    inputs:
        prefix:      string; filename prefix, e.g. "sh_pop_ex"
        file_range:  integer; number of files in consideration, e.g "sh_pop_ex0" to "sh_pop_ex25" is file_range=26
    outputs:
        GS_growth_#:        a number of "GS_growth_" files are generated (same number as specified in file_range),
                            which correspond to the ground-state growth for each state specified in the PYXAID FSSH input
        ex_state_decay_#:   Same as GS_growth. Creates files which correspond to the excited state decays WRT time for each state specified
        NOTE: while GS_growth contains only two columns, one for time and one for GS population, the
              ex_state_decay files contain however many other excited states available.
    '''
    
    GS_growth = []
    ex_state_decay = []
    
    for index in range(0, file_range+1):
        filename = prefix + str(index)
        with open(filename) as f:
            for line in f: # specify index range now for the desired values
                split_line = line.split()
                
                