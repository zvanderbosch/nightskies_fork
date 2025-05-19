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
from astropy import constants

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

def angular_separation(lon1, lat1, lon2, lat2):
    '''
    Compute great circle distance in kilometers between a
    list of coordinates (lon1, lat1) and a single reference
    coordinate (lon2, lat2) using the Vincenty formula:
    
    https://en.wikipedia.org/wiki/Great-circle_distance

    Parameters
    ----------
    lon1: array
        An array of longitude values [radians]
    lat1: array
        An array of latitude values [radians]
    lon2: float
        The reference longitude [radians]
    lat2: float
        The reference latitude [radians]

    Returns
    -------
    dist: array
        Great circle distance [km]
    '''

    # Calculate portions of the Vincenty equation
    deltaRA = abs(lon1 - lon2)
    t1 = n.cos(lat2) * n.sin(deltaRA)
    t2 = n.cos(lat1) * n.sin(lat2)
    t3 = n.sin(lat1) * n.cos(lat2) * n.cos(deltaRA)
    t4 = n.sin(lat1) * n.sin(lat2)
    t5 = n.cos(lat1) * n.cos(lat2) * n.cos(deltaRA)

    # Vincenty equation
    sep = n.arctan2(n.sqrt(t1**2 + (t2-t3)**2), t4+t5)

    # Convert angular separation into distance along Earth's surface
    earthRadius = constants.R_earth.value / 1000
    dist = earthRadius * sep

    return dist


def bearing_angle(lon1, lat1, lon2, lat2):
    '''
    Calculate the bearing angle of (lat2, lon2) with respect to (lat1, lon1). 
    The bearing angle ranges from 0 to 360 degrees, with zero at due north 
    and increasing clockwise. Both the inputs and outputs are in degrees. 

    Parameters
    ----------
    lon1: array
        An array of longitude values [radians]
    lat1: array
        An array of latitude values [radians]
    lon2: float
        The reference longitude [radians]
    lat2: float
        The reference latitude [radians]

    Returns
    -------
    bearing: array
        Bearing angles [deg]

    '''

    # Calculate individual terms
    x = n.cos(lat1)*n.sin(lat2) - n.sin(lat1)*n.cos(lat2)*n.cos(lon1-lon2)
    y = n.sin(lon1-lon2)*n.cos(lat2)

    # Calculate bearing angle
    bearing = n.rad2deg(n.arctan(y/x))

    # Wrap at 180 degrees
    for i in range(len(x)):
        if x[i] < 0: 
            bearing[i] += 180
        bearing[i] = -bearing[i] % 360.
    
    return bearing

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

    # Load in the Places21k spreadsheet with 2010 census data
    placesFile = f"{filepath.scripts}ACP/spreadsheets/Places21k.xlsx"
    places = pd.read_excel(placesFile)
    print(places.head())

    # Calculate great circle distances
    dist = angular_separation(
        n.deg2rad(places.LONGITUDE.values),
        n.deg2rad(places.LATITUDE.values),
        n.deg2rad(lon),
        n.deg2rad(lat)
    )

    # Shorten list of places to those nearby to site
    placesNearby = places.loc[dist < 451].reset_index(drop=True)
    placesNearby['DISTANCE'] = dist[dist < 451]

    # Calculate bearing angles
    bearingAngles = bearing_angle(
        n.deg2rad(placesNearby.LONGITUDE.values),
        n.deg2rad(placesNearby.LATITUDE.values),
        n.deg2rad(lon),
        n.deg2rad(lat)
    )
    placesNearby['BEARING'] = bearingAngles

    # Calculate Walker's law values
    walkersLaw = 0.1 * placesNearby['2010 POPULATION'] * placesNearby['DISTANCE']**(-2.5)
    placesNearby['Walkers Law'] = walkersLaw

    # Shorten list of places to those with Walkers Law > 0.001
    placesHighImpact = placesNearby.loc[walkersLaw > 0.001].reset_index(drop=True)
    print(len(placesHighImpact))



    print(placesHighImpact.head())