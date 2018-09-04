#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  9 15:34:04 2018

@author: oah
"""

def bent_coords(N1, N2, Ag):
    '''
    Ag-N1-N2
    '''
    r1 = abs(Ag-N1)
    r2 = abs(Ag-N2)
    
    #y1 = r1*math.sin(math.pi/6)
    #y2 = r2*math.sin(math.pi/6)
    
    y1 = r1*.5
    y2 = r2*.5
    
    #z1 = y1/math.tan(math.pi/6)
    #z2 = y2/math.tan(math.pi/6)
    
    z1 = r1*0.866025403784439
    z2 = r2*0.866025403784439
    
    print('7    0    0.0  ', str(y1), '    ', str(z1 + Ag))
    print('7    0    0.0  ', str(y2), '    ', str(z2 + Ag))
    
