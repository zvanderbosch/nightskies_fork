#-----------------------------------------------------------------------------#
#albedomodel.py
#
#NPS Night Skies Program
#
#Last updated: 2025/07/24
#
#This script extracts the site-specific albedo value for an observing
#location from a pre-generated geo-referenced albedo model of the
#continental U.S.
#
#Note: 
#
#Input:
#   (1) ws_albedo
#           Raster dataset providing albedo values for continental U.S.
#           (filepath.rasters)
#   (2) wsa_20020610
#           Secondary albedo raster, only used for its coordinate system
#           (filepath.rasters)
#
#Output:
#   (1) siteAlbedo
#           Site-specific albedo value
#
#History:
#	Zach Vanderbosch -- Created script (translated from secondbatchv4.py)
#
#-----------------------------------------------------------------------------#

from glob import glob
from astropy.io import fits

import sys
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
scriptName = 'albedomodel.py'
PREFIX = f'{pc.GREEN}{scriptName:19s}{pc.END}: '

#------------------------------------------------------------------------------#
#-------------------            Define Functions            -------------------#
#------------------------------------------------------------------------------#

def clear_memory(objectList):
    '''
    Function for clearing variables from memory

    Parameters:
    -----------
    objectList: list
        List of objects to delete from memory
    '''
    for obj in objectList:
        if obj:
            del obj


#------------------------------------------------------------------------------#
#-------------------              Main Program              -------------------#
#------------------------------------------------------------------------------#

def calculate_albedo_model(dnight):
    '''
    Main program for determining the site-specific albedo value

    Parameters:
    -----------
    dnight: string
        Name of data night to process

    Returns:
    --------
    siteAlbedo: float
        Site-specific albedo
    '''

    # Filter paths
    F = {'V':'', 'B':'B/'}

    # Print status update
    print(f'{PREFIX}Processing {dnight}...')

    # Set ArcGIS working directories
    scratchsetp = f"{filepath.rasters}scratch_metrics/"
    arcpy.env.workspace = scratchsetp
    arcpy.env.scratchWorkspace = scratchsetp

    # Load in the albedo raster
    albedoRaster = arcpy.sa.Raster(f"{filepath.rasters}ws_albedo")

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

    # Create a point shapefile from longitude and latitude
    pointShapefile = "location1.shp"
    point = arcpy.Point()
    point.X = lon
    point.Y = lat
    pointGeoms = [arcpy.PointGeometry(point)]
    arcpy.management.CopyFeatures(pointGeoms, pointShapefile)

    # Define coordinate system of point shapefile using albedo raster
    dsc = arcpy.Describe(f"{filepath.rasters}wsa_20020610")
    arcpy.management.DefineProjection(
        pointShapefile, dsc.spatialReference
    )

    # Get albedo value from raster for point location
    arcpy.sa.ExtractMultiValuesToPoints(
        pointShapefile, [[albedoRaster, "albedo"]], "NONE"
    )
    rows = arcpy.SearchCursor(pointShapefile)
    row = rows.next()
    siteAlbedo = row.getValue("albedo")
    clear_memory([row,rows])
    print(f'{PREFIX}Site Albedo = {siteAlbedo:.3f}')
    
    return siteAlbedo