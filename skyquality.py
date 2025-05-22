#-----------------------------------------------------------------------------#
#skyquality.py
#
#NPS Night Skies Program
#
#Last updated: 2025/05/22
#
#This script computes sky quality index (SQI) and synthetic
#sky quality meter (SQM) values.
#
#Note: 
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

from glob import glob
from astropy.io import fits

import os
import sys
import pandas as pd

# Local Source
import filepath
import printcolors as pc

# Define print staus prefix
scriptName = 'skyquality.py'
PREFIX = f'{pc.GREEN}{scriptName:19s}{pc.END}: '

#------------------------------------------------------------------------------#
#-------------------            Define Functions            -------------------#
#------------------------------------------------------------------------------#

def get_site_info(imageFile):
    '''
    Function to grab observing site and time info
    '''

    # Load image header
    H = fits.getheader(imageFile)

    # Get site and time info
    lon = H['LONGITUD']
    lat = H['LATITUDE']
    isot = H['DATE-OBS']
    date = isot.split("T")[0]
    time = isot.split("T")[1]

    return lon, lat, date, time


def calc_SQI(gridPath):

    # Load in SQI table
    sqiTableFile = f"{gridPath}sqitbl.dbf"
    sqiTable = pd.read_csv(sqiTableFile)
    print(sqiTable)


#------------------------------------------------------------------------------#
#-------------------              Main Program              -------------------#
#------------------------------------------------------------------------------#

def calculate_sky_quality(dnight,sets):
    '''
    Main program for computing SQI and SQM metrics
    '''

    # Filter paths
    F = {'V':'', 'B':'B/'}


    # Loop through each data set
    for s in sets:

        # Set path for grid datasets
        setnum = int(s[0])
        calsetp = f"{filepath.calibdata}{dnight}/S_{setnum:02d}/{F['V']}"
        gridsetp = f"{filepath.griddata}{dnight}/S_{setnum:02d}/{F['V']}"

        # Get Zenith RA and Dec at dataset midpoint in time
        midpointImage = f"{calsetp}ib022.fit"
        lon,lat,date,time = get_site_info(midpointImage)

        # Calculate sky quality index
        calc_SQI(gridsetp)

