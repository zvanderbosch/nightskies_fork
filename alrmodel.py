#-----------------------------------------------------------------------------#
#alrmodel.py
#
#NPS Night Skies Program
#
#Last updated: 2025/05/19
#
#This script computes All-sky Light pollution Ratio (ALR) model
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
import stat
import arcpy

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
scriptName = 'alrmodel.py'
PREFIX = f'{pc.GREEN}{scriptName:19s}{pc.END}: '

#------------------------------------------------------------------------------#
#-------------------            Define Functions            -------------------#
#------------------------------------------------------------------------------#

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


def clear_memory(objectList):
    '''
    Function for clearing variables from memory
    '''
    for obj in objectList:
        if obj:
            del obj


#------------------------------------------------------------------------------#
#-------------------              Main Program              -------------------#
#------------------------------------------------------------------------------#

def calculate_alr_model(dnight,sets,filter):
    '''
    Main program for computing the numnber/fraction of visible stars
    '''

    # Filter paths
    F = {'V':'', 'B':'B/'}

    # Set ArcGIS working directories
    scratchsetp = f"{filepath.rasters}scratch_metrics/"
    arcpy.env.workspace = scratchsetp
    arcpy.env.scratchWorkspace = scratchsetp

    # Create or clear out working directory
    if os.path.exists(scratchsetp):
        clear_scratch(scratchsetp)
    else:
        os.makedirs(scratchsetp)

    # Load in the ALR raster
    alrRaster = arcpy.sa.Raster(f"{filepath.rasters}alrmodel")

    # Get site longitude and latitude
    imsetp = f"{filepath.calibdata}{dnight}/S_*/{F[filter]}"
    imageFiles = glob(f"{imsetp}ib???.fit")
    if len(imageFiles) == 0:
        print(f'Could not find FITS files at {imsetp}')
        sys.exit(1)
    else:
        with fits.open(imageFiles[0]) as hdul:
            H = hdul[0].header
            lon = H['LONGITUD']
            lat = H['LATITUDE']

    # Create a point shapefile from longitude and latitude
    pointShapefile = "location.shp"
    point = arcpy.Point()
    point.X = lon
    point.Y = lat
    pointGeoms = [arcpy.PointGeometry(point)]
    arcpy.management.CopyFeatures(pointGeoms, pointShapefile)

    # Define coordinate system of point shapefile using alrmodel raster
    dsc = arcpy.Describe(f"{filepath.rasters}alrmodel")
    arcpy.management.DefineProjection(
        pointShapefile, dsc.spatialReference
    )

    # Get ALR value from raster for point location
    arcpy.sa.ExtractMultiValuesToPoints(
        pointShapefile, [[alrRaster, "alr"]], "NONE"
    )
    rows = arcpy.SearchCursor(pointShapefile)
    row = rows.next()
    siteALR = row.getValue("alr")
    clear_memory([row,rows])
    print(f'{PREFIX}Site ALR = {siteALR}')