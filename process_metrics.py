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


def update_progressbar(x,y,t=0):
    '''
    update the progress bar plot
    '''
    if t == 0: 
        #gray out for the in-progress status 
        barax.pcolor([x+4,x+5],[y+5,y+6],[[4],],cmap='gray',vmin=0,vmax=5)
    else:               
        #update for the completed status
        t/=60  # convert to minutes

        #set number formatting
        if t < 10: texty = '%.1f' %t
        else: texty = '%.0f' %t

        #set text color
        if t < 6: c = 'k'
        else: c = 'w'

        #record the time [min] in a master array
        Z[y+5,x] = t

        #set text and background color in figure
        barax.pcolor([x+4,x+5],[y+5,y+6],[[t],],cmap='Greens',vmin=0,vmax=10)
        barax.text(x+4.5, y+5.5, texty, color=c, horizontalalignment='center',
                   verticalalignment='center', size='medium')
        
    #draw the new data and run the GUI's event loop
    plt.pause(0.05)


def process_skyglow(*args):
    '''Calculate luminance/illuminance metrics for artificial skyglow'''
    i = args[-1]
    update_progressbar(0,i)
    t1 = time.time()
    import skyglow as SG
    for filter in args[2]:
        SG.calculate_statistics(args[0],args[1],filter)
    t2 = time.time()
    update_progressbar(0,i,t2-t1)
    print(f'{PREFIX}Processing Time (skyglow): {t2-t1:.2f} seconds')


def process_illumall(*args):
    '''Calculate luminance/illuminance metrics for all light sources'''
    i = args[-1]
    update_progressbar(1,i)
    t1 = time.time()
    import illumall as IA
    for filter in args[2]:
        IA.calculate_statistics(args[0],args[1],filter)
    t2 = time.time()
    update_progressbar(1,i,t2-t1)
    print(f'{PREFIX}Processing Time (illumall): {t2-t1:.2f} seconds')


def process_starsvis(*args):
    '''Calculate visible stars'''
    i = args[-1]
    update_progressbar(2,i)
    t1 = time.time()
    import starsvis as SV
    for filter in args[2]:
        SV.calculate_stars_visible(args[0],args[1],filter)
    t2 = time.time()
    update_progressbar(2,i,t2-t1)
    print(f'{PREFIX}Processing Time (starsvis): {t2-t1:.2f} seconds')


def process_alrmodel(*args):
    '''Calculate All-Sky Light Pollution Ratio (ALR) model'''
    i = args[-1]
    update_progressbar(3,i)
    t1 = time.time()
    import alrmodel as AM
    AM.calculate_alr_model(*args[:-1])
    t2 = time.time()
    update_progressbar(3,i,t2-t1)
    print(f'{PREFIX}Processing Time (alrmodel): {t2-t1:.2f} seconds')


def process_albedomodel(*args):
    '''Calculate albedo model'''
    i = args[-1]
    update_progressbar(4,i)
    t1 = time.time()
    import albedomodel as BM
    BM.calculate_albedo_model(*args[:-1])
    t2 = time.time()
    update_progressbar(4,i,t2-t1)
    print(f'{PREFIX}Processing Time (albedomodel): {t2-t1:.2f} seconds')


def process_places(*args):
    '''Calculate distance & Walker's Law for places'''
    i = args[-1]
    update_progressbar(5,i)
    t1 = time.time()
    import places as PL
    PL.calculate_places(*args[:-1])
    t2 = time.time()
    update_progressbar(5,i,t2-t1)
    print(f'{PREFIX}Processing Time (places): {t2-t1:.2f} seconds')


def process_skyquality(*args):
    '''Calculate SQI and SQM sky quality metrics'''
    i = args[-1]
    update_progressbar(6,i)
    t1 = time.time()
    import skyquality as SQ
    SQ.calculate_sky_quality(*args[:-1])
    t2 = time.time()
    update_progressbar(6,i,t2-t1)
    print(f'{PREFIX}Processing Time (places): {t2-t1:.2f} seconds')


def process_drawmaps(*args):
    '''Generate panoramic graphics'''
    i = args[-1]
    update_progressbar(7,i)
    t1 = time.time()
    import drawmaps as DM
    DM.generate_graphics(*args[:-1])
    t2 = time.time()
    update_progressbar(7,i,t2-t1)
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

    # Read in the processing dataset list and the calibration file names 
    filelist = n.loadtxt(filepath.processlist+'filelist.txt', dtype=str, ndmin=2)
    Dataset, V_band, B_band, _, _, _, processor, centralAz, location = filelist.T
    
    # Determine the number of data sets collected in each night 
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

    # Plot the progress bar template
    barfig, barax = progressbars.bar_metrics(Dataset, nsets)
    print('You have 5 seconds to adjust the position of the progress bar window')
    plt.pause(5) #users have 5 seconds to adjust the figure position

    # Progress bar array (to be filled with processing time)
    Z = n.empty((5+len(filelist),14))*n.nan
    
    
    #------------ Main data processing code --------------------------------------#

    # Looping through multiple data nights
    for i in range(len(filelist)):
        
        # Generate inputs for each processing step
        Filter = []
        if V_band[i] == 'Yes': 
            Filter.append('V')
        if B_band[i] == 'Yes': 
            Filter.append('B')
        
        sets = dnight_sets[Dataset[i]]
        K0 = (Dataset[i],i)
        K1 = (Dataset[i],sets,i)
        K2 = (Dataset[i],sets,Filter,i)
        K3 = (Dataset[i],sets,processor[0],int(centralAz[0]),location[0],i)

        # Status update
        print(
            f'{PREFIX}Processing the {pc.BOLD}{pc.CYAN}'
            f'{Dataset[i]}{pc.END}{pc.END} dataset'
        )

        # Create multiprocessing objects for each step
        p1 = Process(target=process_skyglow,args=K2)
        p2 = Process(target=process_illumall,args=K2)
        p3 = Process(target=process_starsvis,args=K2)
        p4 = Process(target=process_alrmodel,args=K0)
        p5 = Process(target=process_albedomodel,args=K0)
        p6 = Process(target=process_places,args=K0)
        p7 = Process(target=process_skyquality,args=K1)
        p8 = Process(target=process_drawmaps,args=K3)

        # Execute each processing step
        # p1.start()  # Anthropogenic skyglow luminance & illuminance
        # p1.join()
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
        p7.start()  # Sky quality metrics
        p7.join()
        # p8.start()  # Draw maps
        # p8.join()

        # Save the timing records for running the script
        n.savetxt(filepath.calibdata+Dataset[i]+'/processtime_metrics.txt', Z, fmt='%4.1f')
        barfig.savefig(filepath.calibdata+Dataset[i]+'/processtime_metrics.png')

    