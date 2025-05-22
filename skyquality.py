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

from astropy.io import fits
from dbfread import DBF

import numpy as n
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


def calc_SQI(gridPath,mask):

    # Set static weighting factors
    weightingNaturalSky = n.array(
        [0,2,3.2,5,7.9,10,10,10,10,10,10,10,10,10,
         10,10,10,10,10,10,10,10,10,10,10,10,10,10],
        dtype=n.float64
    )
    weightingMilkyWay = n.array(
        [0,1.3,2.1,3.3,5.2,6.65,7.05,7.3,7.7,8.35,9.3,10,
         10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10],
        dtype=n.float64
    )
    weightingScotopicVision = n.array(
        [0,0,0,0,0,0,0,0,0,0,0.55,1.1,1.65,2.2,2.75,3.3,
         3.85,4.4,4.95,5.5,6.05,6.6,7.15,7.7,8.25,8.8,9.35,10],
        dtype=n.float64
    )
    weightingStarsVisible = n.array(
        [0,0.1,0.25,0.45,0.87,1.25,1.7,2.5,3.75,4.85,5.8,6.75,7.35,
         7.85,8.2,8.5,8.8,9,9.2,9.35,9.5,9.6,9.7,9.8,9.87,9.92,9.97,10], 
        dtype=n.float64
    )

    # Load in SQI table corresponding to mask
    if mask == 'horizon':
        sqiTableFile = f"{gridPath}sqitbl.dbf"
    elif mask == 'ZA80':
        sqiTableFile = f"{gridPath}sqitbl80.dbf"
    elif mask == 'ZA70':
        sqiTableFile = f"{gridPath}sqitbl70.dbf"
    sqiTableDBF = DBF(sqiTableFile)
    sqiTable = pd.DataFrame(iter(sqiTableDBF))

    # Get histogram values and order from darkest to brightest bins
    histValues = sqiTable['Value_0'].values[::-1]
    sumValues = sqiTable['Value_0'].sum()

    # Calculate frequency of values in each bin
    frequency = 100 * histValues / sumValues

    # Calculate component linear index values
    indexNaturalSky = 0.25 * sum(weightingNaturalSky * frequency) / 10
    indexMilkyWay = 0.25 * sum(weightingMilkyWay * frequency) / 10
    indexStarsVisible = 0.25 * sum(weightingStarsVisible * frequency) / 10

    # Calculate component log index for scotopic vision
    linearScotopicSum = sum(weightingScotopicVision * frequency)
    if linearScotopicSum <= 0:
        logScotopicSum = 0.0
    else:
        logScotopicSum = n.log10(linearScotopicSum + 1.0)
    logFactor = 100. / n.log10(1000.)
    indexScotopicVision = 0.25 * logFactor * logScotopicSum

    # Calculate final hybrid SQI
    sqi = 100 - indexNaturalSky - indexMilkyWay - indexScotopicVision - indexStarsVisible

    return sqi    


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
        sqiAllsky = calc_SQI(gridsetp, 'horizon')
        sqiZ80 = calc_SQI(gridsetp, 'ZA80')
        sqiZ70 = calc_SQI(gridsetp, 'ZA70')

