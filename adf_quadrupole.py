#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 13:31:59 2018

@author: oah
"""

import numpy as np
import matplotlib.pyplot as plt

for i in ['dip', 'quad', 'tot']:
    for j in ['gauss', 'stick']:
        data = np.genfromtxt('abs_' + j + '_' + i +'.data')
        plt.plot(data[:,0],data[:,1])
    plot_name = i
    plt.title('%s' %plot_name)
    plt.show()
        