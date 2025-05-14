#-----------------------------------------------------------------------------#
#skyglow.py
#
#NPS Night Skies Program
#
#Last updated: 2025/05/14
#
#This script calculates horizontal and vertical illuminance values
#from the anthropogenic light all-aky mosaics.
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
from glob import glob, iglob

import os
import stat
import arcpy
import shutil
import matplotlib.pyplot as plt
import numpy as n

# Local Source
import filepath
import printcolors as pc

# Set arcpy environment variables
arcpy.env.overwriteOutput = True
arcpy.env.rasterStatistics = "NONE"
arcpy.env.pyramid = "NONE"
arcpy.env.compression = "NONE"
arcpy.CheckOutExtension("Spatial")

# Define print staus prefix
PREFIX = f'{pc.GREEN}skyglow.py      {pc.END}: '

#-----------------------------------------------------------------------------#

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


def calculate_illuminance(dnight,sets):
    '''
    Calculate illumninace of anthropogenic skyglow mosaic
    '''

    # Define some paths of interest
    calsetp = f"{filepath.calibdata}{dnight}"
    gridsetp = f"{filepath.griddata}{dnight}"
    scratchsetp = f"{filepath.rasters}scratch_metrics/"

    # Set ArcGIS working directories
    arcpy.env.workspace = scratchsetp
    arcpy.env.scratchWorkspace = scratchsetp

    # Create or clear out scratch directory
    if os.path.exists(scratchsetp):
        clear_scratch(scratchsetp)
    else:
        os.makedirs(scratchsetp)

    # Define paths files
    zoneFile = f'{filepath.rasters}shapefiles/allbands.shp'
    maskRaster = f"{gridsetp}mask/maskd.tif"
