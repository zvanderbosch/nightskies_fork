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
import multiprocessing
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
    args[-1].put(t2-t1)
    print(f'{PREFIX}Processing Time (skyglow): {t2-t1:.2f} seconds')


def process_illumall(*args):
    '''Calculate luminance/illuminance metrics for all light sources'''
    t1 = time.time()
    import illumall as IA
    for filter in args[2]:
        IA.calculate_statistics(args[0],args[1],filter)
    t2 = time.time()
    args[-1].put(t2-t1)
    print(f'{PREFIX}Processing Time (illumall): {t2-t1:.2f} seconds')


def process_starsvis(*args):
    '''Calculate visible stars'''
    t1 = time.time()
    import starsvis as SV
    for filter in args[2]:
        SV.calculate_stars_visible(args[0],args[1],filter)
    t2 = time.time()
    args[-1].put(t2-t1)
    print(f'{PREFIX}Processing Time (starsvis): {t2-t1:.2f} seconds')


def process_alrmodel(*args):
    '''Calculate All-Sky Light Pollution Ratio (ALR) model'''
    t1 = time.time()
    import alrmodel as AM
    AM.calculate_alr_model(*args[:-1])
    t2 = time.time()
    args[-1].put(t2-t1)
    print(f'{PREFIX}Processing Time (alrmodel): {t2-t1:.2f} seconds')


def process_albedomodel(*args):
    '''Calculate albedo model'''
    t1 = time.time()
    import albedomodel as BM
    BM.calculate_albedo_model(*args[:-1])
    t2 = time.time()
    args[-1].put(t2-t1)
    print(f'{PREFIX}Processing Time (albedomodel): {t2-t1:.2f} seconds')


def process_places(*args):
    '''Calculate distance & Walker's Law for places'''
    t1 = time.time()
    import places as PL
    PL.calculate_places(*args[:-1])
    t2 = time.time()
    args[-1].put(t2-t1)
    print(f'{PREFIX}Processing Time (places): {t2-t1:.2f} seconds')


def process_drawmaps(*args):
    '''Generate panoramic graphics'''
    t1 = time.time()
    import drawmaps as DM
    DM.generate_graphics(*args[:-1])
    t2 = time.time()
    args[-1].put(t2-t1)
    print(f'{PREFIX}Processing Time (drawmaps): {t2-t1:.2f} seconds')





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
    Dataset, V_band, B_band, _, _, _, processor, centralAz, location = filelist.T
    
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
        K0 = (Dataset[i],)  
        K1 = (Dataset[i],sets,Filter)  
        K2 = (Dataset[i],sets,processor[0],int(centralAz[0]),location[0])  

        # Status update
        print(
            f'{PREFIX}Processing the {pc.BOLD}{pc.CYAN}'
            f'{Dataset[i]}{pc.END}{pc.END} dataset'
        )

        # Create multiprocessing objects for each step
        q1=Queue(); Q1=(q1,); p1 = Process(target=process_skyglow,args=K1+Q1)
        q2=Queue(); Q2=(q2,); p2 = Process(target=process_illumall,args=K1+Q2)
        q3=Queue(); Q3=(q3,); p3 = Process(target=process_starsvis,args=K1+Q3)
        q4=Queue(); Q4=(q4,); p4 = Process(target=process_alrmodel,args=K0+Q4)
        q5=Queue(); Q5=(q5,); p5 = Process(target=process_albedomodel,args=K0+Q5)
        q6=Queue(); Q6=(q6,); p6 = Process(target=process_places,args=K0+Q6)
        # q7=Queue(); Q7=(q7,); p7 = Process(target=process_sqi,args=K0+Q7)
        q8=Queue(); Q8=(q8,); p8 = Process(target=process_drawmaps,args=K2+Q8)

        # Execute each processing step
        p1.start()  # Anthropogenic skyglow luminance & illuminance
        p1.join()
        # p2.start()  # All sources skyglow luminance & illuminance
        # p2.join()
        # p3.start()  # Number/fraction of visible stars
        # p3.join()
        # p4.start()  # All-sky Light Pollution Ratio (ALR) model
        # p4.join()
        # p5.start()  # Albedo model
        # p5.join()
        # p6.start()  # Places
        # p6.join()
        # p7.start()  # SQI
        # p7.join()
        # p8.start()  # Draw maps
        # p8.join()

    