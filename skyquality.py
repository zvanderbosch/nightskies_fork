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
from photutils.aperture import aperture_photometry
from photutils.aperture import CircularAperture
from photutils.aperture import CircularAnnulus
from photutils.aperture import ApertureStats
from astropy.stats import SigmaClip

import os
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


def measure_skybrightness(imgPath):
    '''
    Function to perform sky brightness aperture
    photometry for SQM calculations.
    '''

    # Load in skybright positions spreadsheet
    allPositions = pd.read_excel(
        f"{filepath.spreadsheets}skybright_positions.xlsx"
    )

    # Iterate over each image
    allPhot = []
    for i in range(45):

        # Get image data
        imgFile = f"{imgPath}ib{i+1:03d}.fit"
        if not os.path.isfile(imgFile):
            continue
        with fits.open(imgFile) as hdul:
            imgData = hdul[0].data

        # Get pixel positions for given image
        imgPositions = allPositions.loc[allPositions.image == i+1]
        xpix = imgPositions.PixelX.values
        ypix = imgPositions.PixelY.values
        xyPositions = [(x,y) for x,y in zip(xpix,ypix)]

        # Perform aperture photometry on image
        apRadius = 20.
        apertures = CircularAperture(xyPositions, r=apRadius)
        photTable = aperture_photometry(imgData, apertures)
        photTable = photTable.to_pandas()

        # Calculate aperture median
        sigclip = SigmaClip(sigma=5.0, maxiters=10)
        aperStats = ApertureStats(imgData, apertures, sigma_clip=sigclip)
        imgPositions.insert(2,'ADU',aperStats.median)

        # Save photometry to list
        allPhot.append(imgPositions)

    # Concatenate photometry into single dataframe
    allPhot = pd.concat(allPhot,ignore_index=True)
    
    return allPhot


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


def calc_SQM(dataNight, setNum, filterName):
    
    # Get zeropoint, extinction coeff, plate scale, & exposure time
    extfile = f"{filepath.calibdata}{dataNight}/extinction_fit_{filterName}.txt"
    extData = n.loadtxt(extfile, ndmin=2)
    zeropoint = extData[setNum-1,2]
    extcoeff = abs(extData[setNum-1,4])
    platescale = extData[setNum-1,8]
    exptime = abs(extData[setNum-1,9])
    psa = 2.5*n.log10((platescale*60)**2) # platescale adjustment

    # Get Zenith RA and Dec at dataset midpoint in time
    imgsetp = f"{filepath.calibdata}{dataNight}/S_{setNum:02d}/"
    midpointImage = f"{imgsetp}ib022.fit"
    lon,lat,date,time = get_site_info(midpointImage)

    # Perform sky brightness measurements
    print(f"{PREFIX}Measuring sky brightness...")
    photometry = measure_skybrightness(imgsetp)

    # Calculate cosine of Zenith angle and apply to ADU
    photometry['cosZ'] = n.cos(n.deg2rad(90.0 - photometry.Pany.values))
    photometry['net_arcsec'] = photometry.cosZ * 3600 * 3600 * 4
    photometry['cosADU'] = photometry.cosZ * photometry.ADU

    # Get sums for zenith angles <= 54 degrees
    arcsecSum = photometry[photometry.Pany >= 36].net_arcsec.sum()
    aduSum = photometry[photometry.Pany >= 36].cosADU.sum()

    # Calculate SQM scale factor
    sqmScaleFactor = 2.5 * n.log10(arcsecSum)

    # Convert platescale adjustment to sq. arcsec per ADU
    psaADU = (4*3600*3600) / (10**(0.4*psa))

    # Calculate synthetic SQM
    mags = zeropoint - 2.5 * n.log10(aduSum * psaADU / exptime)
    sqm = sqmScaleFactor + mags

    return sqm

    
#------------------------------------------------------------------------------#
#-------------------              Main Program              -------------------#
#------------------------------------------------------------------------------#

def calculate_sky_quality(dnight,sets,filter):
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

        # Calculate sky quality indices for each mask
        # sqiAllsky = calc_SQI(gridsetp, 'horizon')
        # sqiZ80 = calc_SQI(gridsetp, 'ZA80')
        # sqiZ70 = calc_SQI(gridsetp, 'ZA70')

        # Calculate SQM
        calc_SQM(dnight, setnum, filter)

