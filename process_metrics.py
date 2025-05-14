#-----------------------------------------------------------------------------#
#process_metrics.py
#
#NPS Night Skies Program
#
#Last updated: 2025/05/14
#
#This script...
#	(1) 
#
#Note: 
#
#
#Input: 
#   (1) 
#
#Output:
#   (1) 
#
#History:
#	Zach Vanderbosch -- Created script
#
#-----------------------------------------------------------------------------#

import os
import time
import warnings
import numpy as n
import matplotlib.pyplot as plt

from datetime import datetime as Dtime
from multiprocessing import Process, Queue

# Local source
import filepath
import progressbars
import printcolors as pc

# Define print status prefix
PREFIX = f'{pc.GREEN}process_metrics.py {pc.END}: '


##########################  Definitions  ######################################




##########################  Main Program  #####################################

if __name__ == '__main__':
    t1 = time.time()
    #-----------------------------------------------------------------------------#    
    print(
        '\n--------------------------------------------------------------\n\n'
        '     NPS NIGHT SKIES PROGRAM LIGHT POLLUTION CALCULATIONS           '
        '\n\n--------------------------------------------------------------\n'
    )
        
    #------------ Read in the processing list and initialize ---------------------#

    #Read in the processing dataset list and the calibration file names 
    filelist = n.loadtxt(filepath.processlist+'filelist.txt', dtype=str, ndmin=2)
    Dataset, V_band, B_band, Flat_V, Flat_B, Curve, Processor = filelist.T
    
    #Check the calibration files exist    
    for i in range(len(filelist)):
        if V_band == 'Yes':
            open(filepath.flats+Flat_V[i])
        if B_band == 'Yes':
            open(filepath.flats+Flat_B[i])
        open(filepath.lincurve+Curve[i]+'.txt')
    
    #Determine the number of data sets collected in each night 
    img_sets = set(['1st','2nd','3rd','4th','5th','6th','7th','8th'])
    dnight_sets = {}
    nsets = []
    for dnight in Dataset:
        n_path = filepath.rawdata + dnight + '/'
        dnight_sets[dnight] = []
        for f in os.listdir(n_path):
            if os.path.isdir(n_path+f) & (f in img_sets):
                dnight_sets[dnight].append(f)
        nsets.append(len(dnight_sets[dnight]))
        
        #Make calibration folders
        if not os.path.exists(filepath.calibdata+dnight):
            os.makedirs(filepath.calibdata+dnight)
    
    
    #------------ Main data processing code --------------------------------------#
    