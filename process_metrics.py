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
scriptName = 'process_metrics.py'
PREFIX = f'{pc.GREEN}{scriptName:19s}{pc.END}: '


######################  Function Definitions  #################################


def process_skyglow(*args):
    '''Calculate luminance/illuminance metrics for artificial skyglow'''
    t1 = time.time()
    import skyglow as SG
    for filter in args[2]:
        SG.calculate_statistics(args[0],args[1],filter)
    t2 = time.time()
    print(f'{PREFIX}Processing Time (skyglow): {t2-t1:.2f} seconds')


def process_illumall(*args):
    '''Calculate luminance/illuminance metrics for all light sources'''
    t1 = time.time()
    import illumall as IA
    for filter in args[2]:
        IA.calculate_statistics(args[0],args[1],filter)
    t2 = time.time()
    print(f'{PREFIX}Processing Time (illumall): {t2-t1:.2f} seconds')

def process_starsvis(*args):
    '''Calculate visible stars'''
    t1 = time.time()
    import starsvis as SV
    for filter in args[2]:
        SV.calculate_stars_visible(args[0],args[1],filter)
    t2 = time.time()
    print(f'{PREFIX}Processing Time (starsvis): {t2-t1:.2f} seconds')





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
    Dataset, V_band, B_band, _, _, _, Processor = filelist.T
    
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

    #Looping through multiple data nights
    for i in range(len(filelist)):
        
        # Generate inputs for each processing step
        Filter = []
        if V_band[i] == 'Yes': 
            Filter.append('V')
        if B_band[i] == 'Yes': 
            Filter.append('B')
        
        sets = dnight_sets[Dataset[i]]
        # K0 = (Dataset[i],sets,Filterset,Curve[i])
        K1 = (Dataset[i],sets,Filter) 
        # K2 = (Dataset[i],sets)  

        # Status update
        print(
            f'{PREFIX}Processing the {pc.BOLD}{pc.CYAN}'
            f'{Dataset[i]}{pc.END}{pc.END} dataset'
        )

        # Execute each processing step
        # process_skyglow(*K1)               # Asthropogenic skyglow luminance & illuminance
        # process_illumall(*K1)              # All sources skyglow luminance & illuminance
        # process_starsvis(*K1)              # Number/fraction of visible stars

    