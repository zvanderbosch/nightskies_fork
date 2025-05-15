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
scriptName = 'skyglow.py'
PREFIX = f'{pc.GREEN}{scriptName:19s}{pc.END}: '

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

def clear_memory(objectList):
    '''
    Function for clearing variables from memory
    '''
    for obj in objectList:
        if obj:
            del obj

def nl_to_mlux(raster):
    '''
    Unit conversion from nano-Lamberts to milli-Lux
    '''
    return (0.000000761544 * raster) / 314.159


def calc_sky_luminance(mosaicDict, zoneRaster, gridPath, results):
    '''
    Function for calculating sky luminance metrics
    '''

    # Iterate over each mosaic
    for i,mosaic in enumerate(mosaicDict.values()):

        # Calculate zonal stats
        outputTable = f"{gridPath}skyhemis{i}.dbf"
        _ = arcpy.sa.ZonalStatisticsAsTable(
            zoneRaster,"VALUE",mosaic,outputTable,"DATA","ALL"
        )

        # Extract stats from output table
        rows = arcpy.SearchCursor(outputTable)
        row = rows.next()
        results[f'skymax{i}'] = row.getValue("MAX")
        results[f'skymin{i}'] = row.getValue("MIN")
        results[f'skyave{i}'] = row.getValue("MEAN")
        clear_memory([row,rows])

    return results


def calc_zonal_sky_luminance(mosaic, zoneRaster, gridPath, results):
    '''
    Function for calculating mean sky luminance by zone
    '''

    # Calculate zonal stats
    outputTable = f"{gridPath}skyzones.dbf"
    _ = arcpy.sa.ZonalStatisticsAsTable(
        zoneRaster,"band",mosaic,outputTable,"DATA","MEAN"
    )

    # Extract stats from output table
    rows = arcpy.SearchCursor(outputTable)
    for i in range(5):
        row = rows.next()
        results[f'zoneAve{i}'] = row.getValue("MEAN")
        results[f'zoneMax{i}'] = row.getValue("COUNT")
    clear_memory([row,rows])

    return results


def calc_luminouos_emittance(mosaicDict, zoneRaster, gridPath, results):
    '''
    Function for calculating sky luminance metrics
    '''

    # Convert from nL to mlux units
    mlux = nl_to_mlux(mosaicDict['allsky'])
    mlux80 = nl_to_mlux(mosaicDict['za80'])
    mlux70 = nl_to_mlux(mosaicDict['za70'])
    mosaicListMlux = [mlux, mlux80, mlux70]

    # Iterate over each mosaic
    for i,mosaic in enumerate(mosaicListMlux):

        # Calculate zonal stats
        outputTable = f"{gridPath}skyhemis{i}.dbf"
        _ = arcpy.sa.ZonalStatisticsAsTable(
            zoneRaster,"VALUE",mosaic,outputTable,"DATA","ALL"
        )

        # Extract stats from output table
        rows = arcpy.SearchCursor(outputTable)
        row = rows.next()
        results[f'hemis{i}'] = row.getValue("MEAN")
        results[f'totalill{i}'] = row.getValue("SUM")
        clear_memory([row,rows])

    return results


def calculate_illuminance(dnight,sets,filter):
    '''
    Calculate illumninace of anthropogenic skyglow mosaic
    '''

    # Filter paths
    F = {'V':'', 'B':'B/'}

    # Define path to scratch workspace
    scratchsetp = f"{filepath.rasters}scratch_metrics/"

    # Set ArcGIS working directories
    arcpy.env.workspace = scratchsetp
    arcpy.env.scratchWorkspace = scratchsetp

    # Create or clear out scratch directory
    if os.path.exists(scratchsetp):
        clear_scratch(scratchsetp)
    else:
        os.makedirs(scratchsetp)

    # Define paths to a few needed files
    zoneFile = f'{filepath.rasters}shapefiles/allbands.shp'
    za80File = f'{filepath.rasters}shapefiles/80zamaskf.shp'
    za70File = f'{filepath.rasters}shapefiles/70zamaskf.shp'
    maskRasterFile = f"{filepath.griddata}{dnight}/mask/maskd.tif"

    # Load in the mask raster
    maskRaster = arcpy.sa.Raster(maskRasterFile)

    # Loop through each data set
    for s in sets:

        # Set path for grid datasets
        setnum = int(s[0])
        gridsetp = f"{filepath.griddata}{dnight}/S_{setnum:02d}/{F[filter]}"

        # Initial status update for data set
        print(f'{PREFIX}Processing {dnight} Set-{setnum} {filter}-band...')

        # Load in anthropogenic skyglow mosaic in nano-Lamberts [nL]
        anthRaster = arcpy.sa.Raster(f"{gridsetp}skyglow/anthlightnl")

        # Clip mosaic to 80 and 70 zenith-angle limits
        print(f'{PREFIX}Clipping skyglow mosaic to 80 and 70 Zenith Angle limits...')
        arcpy.management.Clip(
            anthRaster,"#",f"anth80{setnum}",za80File,"#","ClippingGeometry"
        )
        arcpy.management.Clip(
            anthRaster,"#",f"anth70{setnum}",za70File,"#","ClippingGeometry"
        )
        anthRaster80 = arcpy.sa.Raster(f"anth80{setnum}")
        anthRaster70 = arcpy.sa.Raster(f"anth70{setnum}")

        # Initialize result and mosaic dicts
        resultDict = {}
        mosaicDict = {
            'allsky': anthRaster,
            'za80': anthRaster80,
            'za70': anthRaster70
        }

        # Calculate sky luminance metrics
        print(f'{PREFIX}Calculating anthropogenic sky luminance...')
        calc_sky_luminance(mosaicDict, maskRaster, gridsetp, resultDict)


        # Calculate mean sky luminance by zone
        print(f'{PREFIX}Calculating mean anthropogenic sky luminance by zone...')
        calc_zonal_sky_luminance(anthRaster, zoneFile, gridsetp, resultDict)


        # Calculate zonal stats for each mosaic
        print(f'{PREFIX}Calculating all-sky anthropogenic luminous emittance...')
        calc_luminouos_emittance(mosaicDict, maskRaster, gridsetp, resultDict)


        # Print out results
        for key in resultDict:
            print(key, resultDict[key])
            
