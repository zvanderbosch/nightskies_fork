#-----------------------------------------------------------------------------#
#places.py
#
#NPS Night Skies Program
#
#Last updated: 2025/05/19
#
#This script calculates great-circle distances and Walker's law values
#to nearby places (cities/towns) using 2010 Census Data.
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
#	Zach Vanderbosch -- Created script (translated from places21g.vbs)
#
#-----------------------------------------------------------------------------#

from glob import glob
from astropy.io import fits

import os
import sys
import stat

# Local Source
import filepath
import printcolors as pc

# Define print staus prefix
scriptName = 'places.py'
PREFIX = f'{pc.GREEN}{scriptName:19s}{pc.END}: '

#------------------------------------------------------------------------------#
#-------------------            Define Functions            -------------------#
#------------------------------------------------------------------------------#



#------------------------------------------------------------------------------#
#-------------------              Main Program              -------------------#
#------------------------------------------------------------------------------#

def calculate_places(dnight):
    '''
    Main program for computing great circle distances and Walker's
    Law values for nearby cities and towns using 2010 Census data.
    '''

    # Filter paths
    F = {'V':'', 'B':'B/'}

    # Get site longitude and latitude
    imsetp = f"{filepath.calibdata}{dnight}/S_*/{F['V']}"
    imageFiles = sorted(glob(f"{imsetp}ib???.fit"))
    if len(imageFiles) == 0:
        print(f'Could not find FITS files at {imsetp}')
        sys.exit(1)
    else:
        with fits.open(imageFiles[0]) as hdul:
            H = hdul[0].header
            lon = H['LONGITUD']
            lat = H['LATITUDE']