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
import sys
import time
import warnings
import numpy as n
import pandas as pd
import matplotlib.pyplot as plt

from datetime import datetime as Dtime
from multiprocessing import Process, Queue

# Add path to ccdmodules
sys.path.append('./ccdmodules')

# Local source
import filepath as filepath
import progressbars as pb
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
    t1 = time.time()
    import ccdmodules.skyglow as SG
    sgResults = []
    for filter in args[2]:
        if filter != "V":
            continue
        sgEntry = SG.calculate_statistics(args[0],args[1],filter)
        sgResults.append(sgEntry)
    sgResults = pd.concat(sgResults,ignore_index=True)
    t2 = time.time()
    args[-1].put([t2-t1,sgResults])
    print(f'{PREFIX}Processing Time (skyglow): {t2-t1:.2f} seconds')


def process_illumall(*args):
    '''Calculate luminance/illuminance metrics for all light sources'''
    t1 = time.time()
    import ccdmodules.illumall as IA
    iaResults = []
    for filter in args[2]:
        if filter != "V":
            continue
        iaEntry = IA.calculate_statistics(args[0],args[1],filter)
        iaResults.append(iaEntry)
    iaResults = pd.concat(iaResults,ignore_index=True)
    t2 = time.time()
    args[-1].put([t2-t1,iaResults])
    print(f'{PREFIX}Processing Time (illumall): {t2-t1:.2f} seconds')


def process_starsvis(*args):
    '''Calculate visible stars'''
    t1 = time.time()
    import ccdmodules.starsvis as SV
    svResults = []
    for filter in args[2]:
        if filter != "V":
            continue
        svEntry = SV.calculate_stars_visible(args[0],args[1],filter)
        svResults.append(svEntry)
    svResults = pd.concat(svResults, ignore_index=True)
    t2 = time.time()
    args[-1].put([t2-t1,svResults])
    print(f'{PREFIX}Processing Time (starsvis): {t2-t1:.2f} seconds')


def process_alrmodel(*args):
    '''Calculate All-Sky Light Pollution Ratio (ALR) model'''
    t1 = time.time()
    import ccdmodules.alrmodel as AM
    alr = AM.calculate_alr_model(*args[:-1])
    t2 = time.time()
    args[-1].put([t2-t1, alr])
    print(f'{PREFIX}Processing Time (alrmodel): {t2-t1:.2f} seconds')


def process_albedomodel(*args):
    '''Calculate albedo model'''
    t1 = time.time()
    import ccdmodules.albedomodel as BM
    albedo = BM.calculate_albedo_model(*args[:-1])
    t2 = time.time()
    args[-1].put([t2-t1, albedo])
    print(f'{PREFIX}Processing Time (albedomodel): {t2-t1:.2f} seconds')


def process_places(*args):
    '''Calculate distance & Walker's Law for places'''
    t1 = time.time()
    import ccdmodules.places as PL
    PL.calculate_places(*args[:-1])
    t2 = time.time()
    args[-1].put(t2-t1)
    print(f'{PREFIX}Processing Time (places): {t2-t1:.2f} seconds')


def process_skyquality(*args):
    '''Calculate SQI and SQM sky quality metrics'''
    t1 = time.time()
    import ccdmodules.skyquality as SQ
    sqResults = []
    for filter in args[2]:
        if filter != "V":
            continue
        sqEntry = SQ.calculate_sky_quality(args[0],args[1],filter,args[3])
        sqResults.append(sqEntry)
    sqResults = pd.concat(sqResults,ignore_index=True)
    t2 = time.time()
    args[-1].put([t2-t1,sqResults])
    print(f'{PREFIX}Processing Time (places): {t2-t1:.2f} seconds')


def process_drawmaps(*args):
    '''Generate panoramic graphics'''
    t1 = time.time()
    import ccdmodules.drawmaps as DM
    DM.generate_graphics(*args[:-1])
    t2 = time.time()
    args[-1].put(t2-t1)
    print(f'{PREFIX}Processing Time (drawmaps): {t2-t1:.2f} seconds')


def process_tables(*args):
    '''Generate summary tables'''
    t1 = time.time()
    import ccdmodules.savetables as ST
    ST.generate_tables(*args[:-1])
    t2 = time.time()
    args[-1].put(t2-t1)
    print(f'{PREFIX}Processing Time (savetables): {t2-t1:.2f} seconds')



##########################  Main Program  #####################################

if __name__ == '__main__':

    t1 = time.time()
    #-----------------------------------------------------------------------------#    
    print(
        '\n--------------------------------------------------------------\n\n'
        '     NPS NIGHT SKIES PROGRAM LIGHT POLLUTION CALCULATIONS           '
        '\n\n--------------------------------------------------------------\n'
    )
    warnings.filterwarnings("ignore",".*GUI is implemented.*")
        
    #------------ Read in the processing list and initialize ---------------------#

    # Read in the processing dataset list and the calibration file names 
    filelist = n.loadtxt(filepath.processlist+'filelist.txt', dtype=str, ndmin=2)
    Dataset, V_band, B_band, _,_,_,_, processor, centralAz, location = filelist.T
    
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
    # barfig, barax = pb.bar_metrics(Dataset, nsets)
    # print('You have 5 seconds to adjust the position of the progress bar window')
    # plt.pause(5) #users have 5 seconds to adjust the figure position

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
        K0 = (Dataset[i],)
        K1 = (Dataset[i],sets)
        K2 = (Dataset[i],sets,Filter)
        K3 = (Dataset[i],sets,processor[0],int(centralAz[0]),location[0])

        # Status update
        print(
            f'{PREFIX}Processing the {pc.BOLD}{pc.CYAN}'
            f'{Dataset[i]}{pc.END}{pc.END} dataset'
        )

        # Anthropogenic skyglow luminance & illuminance
        q0=Queue(); args=(Dataset[i],sets,Filter,q0)
        p0 = Process(target=process_skyglow,args=args)
        p0.start(); #update_progressbar(0,i)
        p0.join() ; #update_progressbar(0,i,q0.get()[0])
        skyglowMetrics = q0.get()[1]
        
        # All sources skyglow luminance & illuminance
        q1=Queue(); args=(Dataset[i],sets,Filter,q1)
        p1 = Process(target=process_illumall,args=args)
        p1.start(); #update_progressbar(1,i)
        p1.join() ; #update_progressbar(1,i,q1.get()[0])
        illumammMetrics = q1.get()[1]

        # # Number/fraction of visible stars
        # q2=Queue(); args=(Dataset[i],sets,Filter,q2)
        # p2 = Process(target=process_starsvis,args=args)
        # p2.start(); #update_progressbar(2,i)
        # p2.join() ; #update_progressbar(2,i,q2.get()[0])
        # numstars = q2.get()[1]

        # # All-sky Light Pollution Ratio (ALR) model
        # q3=Queue(); args=(Dataset[i],q3)
        # p3 = Process(target=process_alrmodel,args=args)
        # p3.start(); #update_progressbar(3,i)
        # p3.join() ; #update_progressbar(3,i,q3.get()[0])
        # siteALR = q3.get()[1]

        # # Calculate Site Albedo
        # q4=Queue(); args=(Dataset[i], q4); 
        # p4 = Process(target=process_albedomodel,args=args)
        # p4.start(); #update_progressbar(4,i)
        # p4.join() ; #update_progressbar(4,i,q4.get()[0])
        # siteAlbedo = q4.get()[1]

        # # Places
        # q5=Queue(); args=(Dataset[i], q5)
        # p5 = Process(target=process_places,args=args)
        # # p5.start(); update_progressbar(5,i)
        # # p5.join() ; update_progressbar(5,i,q5.get())

        # # Sky quality metrics
        # q6=Queue(); args=(Dataset[i],sets,Filter,siteAlbedo,q6)
        # p6 = Process(target=process_skyquality,args=args)
        # p6.start(); #update_progressbar(6,i)
        # p6.join() ; #update_progressbar(6,i,q6.get()[0])
        # sqMetrics = q6.get()[1]

        # # Draw maps
        # q7=Queue(); args=(Dataset[i],sets,processor[0],int(centralAz[0]),location[0],q7)
        # p7 = Process(target=process_drawmaps,args=args)
        # # p7.start(); update_progressbar(7,i)
        # # p7.join() ; update_progressbar(7,i,q7.get())

        # # Combine light pollution metrics in a single parameter
        # metricResults = {
        #     'skyglow': skyglowMetrics,
        #     'starsvis': numstars,
        #     'alr': siteALR,
        #     'albedo': siteAlbedo,
        #     'skyquality': sqMetrics
        # }

        # # Execute save tables process
        # q8=Queue()
        # tableArgs = (
        #     Dataset[i], 
        #     sets, 
        #     processor[0],
        #     int(centralAz[0]),
        #     location[0], 
        #     metricResults,
        #     q8
        # )
        # p8 = Process(target=process_tables,args=tableArgs)
        # p8.start(); #update_progressbar(8,i)
        # p8.join() ; #update_progressbar(8,i,q8.get())

        # # Save the timing records for running the script
        # # n.savetxt(filepath.calibdata+Dataset[i]+'/processtime_metrics.txt', Z, fmt='%4.1f')
        # # barfig.savefig(filepath.calibdata+Dataset[i]+'/processtime_metrics.png')

    