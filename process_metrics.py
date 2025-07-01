#-----------------------------------------------------------------------------#
#process_metrics.py
#
#NPS Night Skies Program
#
#Last updated: 2025/06/24
#
#This script executes the final steps of the NSNSD CCD data processing
#pipeline, analyzing the observed sky brightness and inferred artificial
#light mosaics and generating a variety of light pollution metrics. The 
#final outputs are a table containing a summary of the metrics and metadata
#associated with data collection and calibration, and graphics displaying
#panoramic images of the observed sky conditions, the natural sky model,
#and the inferred contribution to sky brightness from artificial light.
#
#
#Note: 
#
#
#Input: 
#   (1) filelist.xlsx
#   (2) Sky brightness mosaic (skybrightnl)
#   (3) Anthropogenic light mosaic (anthlightnl)
#
#Output:
#   (1) List of nearby cities ranked by Walker's Law values (calibdata)
#       a. cities.xlsx
#   (2) Vertical & Horizontal Illuminance tables & plots (calibdata)
#       a. illuminance_horizon.png
#       b. illuminance_za80.png
#       c. illuminance_za70.png
#       d. vert.xlsx
#   (3) Panoramic graphics (graphics)
#       a. <DATANIGHT>_fullres_<DATASET>_HA<CENTRAL-AZIMUTH>.jpg
#       b. <DATANIGHT>_skybright_<DATASET>_HA<CENTRAL-AZIMUTH>.jpg
#       c. <DATANIGHT>_natsky_<DATASET>_HA<CENTRAL-AZIMUTH>.jpg
#       d. <DATANIGHT>_artificial_<DATASET>_HA<CENTRAL-AZIMUTH>.jpg
#   (4) Sky Brightness and Light Pollution Metrics (tables)
#       a. <DATANIGHT>.xlsx
#   (5) Processing summaries
#       a. processing_metrics.png
#       b. processing_metrics.txt
#
#
#
#History:
#	Zach Vanderbosch -- Created script (translated from secondbatchv4.py)
#
#-----------------------------------------------------------------------------#

import os
import sys
import time
import stat
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
    print(f'{PREFIX}Processing Time (skyquality): {t2-t1:.2f} seconds')


def process_drawgraphics(*args):
    '''Generate panoramic graphics'''
    t1 = time.time()
    import drawgraphics as DG
    DG.generate_graphics(*args[:-1])
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


def clear_scratch(scratch_dir):
    '''
    Function to clear out all files and folders from
    the scratch directory.
    '''
    for root, dirs, files in os.walk(scratch_dir, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.chmod(os.path.join(root, name), stat.S_IWRITE)
            os.rmdir(os.path.join(root, name))


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

    # Read in processing dataset list and skip rows where Process = No
    filelist = pd.read_excel(f"{filepath.processlist}filelist.xlsx")
    filelist = filelist.loc[filelist['Process'] == 'Yes'].reset_index(drop=True)

    # Get metadata for each dataset to be processed
    Dataset = filelist['Dataset'].values
    V_band = filelist['V_band'].values
    B_band = filelist['B_band'].values
    processor = filelist['Processor'].values
    centralAz = filelist['Central_AZ'].values
    location = filelist['Location'].values
    
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
        
        #Alert if calibdata directory not found
        calsetp = f"{filepath.calibdata}{dnight}"
        if not os.path.exists(calsetp):
            print(f'ERROR: Calibdata directory not found at {calsetp}')
            sys.exit(1)

    # Plot the progress bar template
    barfig, barax = pb.bar_metrics(Dataset, nsets)
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
        
        # Get available datasets for the given night
        sets = dnight_sets[Dataset[i]]

        # Status update
        print(
            f'{PREFIX}Processing the {pc.BOLD}{pc.CYAN}'
            f'{Dataset[i]}{pc.END}{pc.END} dataset'
        )

        # Create or clear out scratch directory
        scratchsetp = f"{filepath.rasters}scratch_metrics/"
        if os.path.exists(scratchsetp):
            clear_scratch(scratchsetp)
        else:
            os.makedirs(scratchsetp)

        # Remove outputs from any prior runs
        if os.path.isfile(f"{filepath.calibdata}{Dataset[i]}/vert.xlsx"):
            os.remove(f"{filepath.calibdata}{Dataset[i]}/vert.xlsx")
        if os.path.isfile(f"{filepath.calibdata}{Dataset[i]}/cities.xlsx"):
            os.remove(f"{filepath.calibdata}{Dataset[i]}/cities.xlsx")
        if os.path.isfile(f"{filepath.tables}{Dataset[i]}.xlsx"):
            os.remove(f"{filepath.tables}{Dataset[i]}.xlsx")

        # Setup first six processes
        q0=Queue(); args=(Dataset[i],sets,Filter,q0)
        p0 = Process(target=process_illumall,args=args)     # All sources skyglow luminance & illuminance
        q1=Queue(); args=(Dataset[i],sets,Filter,q1)
        p1 = Process(target=process_skyglow,args=args)      # Anthropogenic skyglow luminance & illuminance
        q2=Queue(); args=(Dataset[i],sets,Filter,q2)
        p2 = Process(target=process_starsvis,args=args)     # Number/fraction of visible stars
        q3=Queue(); args=(Dataset[i],q3)
        p3 = Process(target=process_alrmodel,args=args)     # All-sky Light Pollution Ratio (ALR) model
        q4=Queue(); args=(Dataset[i], q4); 
        p4 = Process(target=process_albedomodel,args=args)  # Calculate Site Albedo
        q5=Queue(); args=(Dataset[i], q5)
        p5 = Process(target=process_places,args=args)       # Places

        # Execute first six processes
        p0.start(); update_progressbar(0,i)
        p1.start(); update_progressbar(1,i)
        p2.start(); update_progressbar(2,i)
        p3.start(); update_progressbar(3,i)
        p4.start(); update_progressbar(4,i)
        p5.start(); update_progressbar(5,i)
        p5.join() ; r5=q5.get(); update_progressbar(5,i,r5)
        p3.join() ; r3=q3.get(); update_progressbar(3,i,r3[0])
        p4.join() ; r4=q4.get(); update_progressbar(4,i,r4[0])
        p2.join() ; r2=q2.get(); update_progressbar(2,i,r2[0])
        p0.join() ; r0=q0.get(); update_progressbar(0,i,r0[0])
        p1.join() ; r1=q1.get(); update_progressbar(1,i,r1[0])

        # Get metrics from returned results
        illumallMetrics = r0[1]
        skyglowMetrics = r1[1]
        numstars = r2[1]
        siteALR = r3[1]
        siteAlbedo = r4[1]

        # Setup the next two processes
        q6=Queue(); args=(Dataset[i],sets,Filter,siteAlbedo,q6)
        p6 = Process(target=process_skyquality,args=args)   # Sky quality metrics
        q7=Queue(); args=(Dataset[i],sets,processor[i],int(centralAz[i]),location[i],q7)
        p7 = Process(target=process_drawgraphics,args=args)     # Draw maps
        
        # Execute next two processes
        p6.start(); update_progressbar(6,i)
        p7.start(); update_progressbar(7,i)
        p6.join() ; r6=q6.get(); update_progressbar(6,i,r6[0])
        p7.join() ; r7=q7.get(); update_progressbar(7,i,r7)

        # Get metrics from returned results
        skyqualityMetrics = r6[1]

        # Clear out scratch directory when finished
        clear_scratch(scratchsetp)

        # Load statistics from natsky_model_params.xlsx
        natskyAllSources = pd.read_excel(
            f"{filepath.calibdata}{Dataset[i]}/natsky_model_params.xlsx",
            sheet_name="Sky_Brightness_All_Sources"
        )
        natskyArtificial = pd.read_excel(
            f"{filepath.calibdata}{Dataset[i]}/natsky_model_params.xlsx",
            sheet_name="Sky_Brightness_Artificial_Only"
        )

        # Combine light pollution metrics in a single parameter
        metricResults = {
            'skyglow': skyglowMetrics,
            'illumall': illumallMetrics,
            'starsvis': numstars,
            'alr': siteALR,
            'albedo': siteAlbedo,
            'skyquality': skyqualityMetrics,
            'natsky_all': natskyAllSources,
            'natsky_art': natskyArtificial
        }

        # Execute save tables process
        q8=Queue()
        tableArgs = (
            Dataset[i], 
            sets, 
            processor[i],
            int(centralAz[i]),
            location[i], 
            metricResults,
            q8
        )
        p8 = Process(target=process_tables,args=tableArgs)
        p8.start(); update_progressbar(8,i)
        p8.join() ; r8 = q8.get()
        update_progressbar(8,i,r8)

        # Save the timing records for running the script
        n.savetxt(filepath.calibdata+Dataset[i]+'/processtime_metrics.txt', Z, fmt='%4.1f')
        barfig.savefig(filepath.calibdata+Dataset[i]+'/processtime_metrics.png')

    