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
import numpy as n
import pandas as pd

# Local Source
import filepath
import printcolors as pc

# Define print staus prefix
scriptName = 'places.py'
PREFIX = f'{pc.GREEN}{scriptName:19s}{pc.END}: '

#------------------------------------------------------------------------------#
#-------------------            Define Functions            -------------------#
#------------------------------------------------------------------------------#

def angular_separation(ra1, de1, ra2, de2):
    '''
    Compute great circle angular separation between a list
    of coordinates (ra1, de1) and a single reference
    coordinate (ra2, de2) using the Vincenty formula:
    
    https://en.wikipedia.org/wiki/Great-circle_distance

    Parameters
    ----------
    ra1: array
        An array of RA values [radians]
    de1: array
        An array of Dec values [radians]
    ra2: float
        The reference RA [radians]
    de2: float
        The reference Dec [radians]

    Returns
    -------
    sep: float
        Angular separation in degrees
    '''

    # Calculate portions of the Vincenty equation
    deltaRA = abs(ra1 - ra2)
    t1 = n.cos(de2) * n.sin(deltaRA)
    t2 = n.cos(de1) * n.sin(de2)
    t3 = n.sin(de1) * n.cos(de2) * n.cos(deltaRA)
    t4 = n.sin(de1) * n.sin(de2)
    t5 = n.cos(de1) * n.cos(de2) * n.cos(deltaRA)

    # Vincenty equation
    sep = n.arctan2(n.sqrt(t1**2 + (t2-t3)**2), t4+t5)
    sep = n.rad2deg(sep)

    return sep

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

    print(lat, lon)

    # Load in the Places21k spreadsheet with 2010 census data
    placesFile = f"{filepath.scripts}ACP/spreadsheets/Places21k.xlsx"
    places = pd.read_excel(placesFile)
    print(places.head())

    # Shorten list of places to those nearby to site
    placesNearby = places[
        (places
    ]