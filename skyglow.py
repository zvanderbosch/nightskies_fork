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

def raster_stats(raster):

    # Convert raster to numpy array
    arr = arcpy.RasterToNumPyArray(
        raster,
        nodata_to_value=n.nan
    )

    # Calculate stats
    statDict = {
        'MEAN': n.nanmean(arr),
        'MIN': n.nanmin(arr),
        'MAX': n.nanmax(arr),
        'SUM': n.nansum(arr)
    }


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
    za70File = f'{filepath.rasters}shapefiles/70zamaskf.shp'
    za80File = f'{filepath.rasters}shapefiles/80zamaskf.shp'
    maskRasterFile = f"{filepath.griddata}{dnight}/mask/maskd.tif"

    # Load in the mask raster
    maskRaster = arcpy.sa.Raster(maskRasterFile)

    # Loop through each data set
    for s in sets:

        # Set path for grid datasets
        setnum = int(s[0])
        gridsetp = f"{filepath.griddata}{dnight}/S_{setnum:02d}/{F[filter]}"

        # Status update
        print(f'{PREFIX}Processing {dnight} Set-{setnum} {filter}-band...')

        # Load in anthropogenic skyglow mosaic in nano-Lamberts [nL]
        anthRaster = arcpy.sa.Raster(f"{gridsetp}skyglow/anthlightnl")

        # Clip mosaic to 70 and 80 zenith-angle limits
        print(f'{PREFIX}Clipping skyglow mosaic to 70 and 80 Zenith Angle limits...')
        arcpy.management.Clip(
            anthRaster,"#",f"anth70{setnum}",za70File,"#","ClippingGeometry"
        )
        arcpy.management.Clip(
            anthRaster,"#",f"anth80{setnum}",za80File,"#","ClippingGeometry"
        )
        anthRaster70 = arcpy.sa.Raster(f"anth70{setnum}")
        anthRaster80 = arcpy.sa.Raster(f"anth80{setnum}")

        # Calculate zonal stats for each mosaic
        print(f'{PREFIX}Calculating anthropogenic sky luminance...')
        resultDict = {}
        mosaicList = [anthRaster,anthRaster70,anthRaster80]
        for i,mosaic in enumerate(mosaicList):

            # Calculate zonal stats
            outputTable = f"{gridsetp}skyhemis{i}.dbf"
            _ = arcpy.sa.ZonalStatisticsAsTable(
                maskRaster,"VALUE",mosaic,outputTable,"DATA","ALL"
            )

            # Extract stats from output table
            rows = arcpy.SearchCursor(outputTable)
            row = rows.next()
            resultDict[f'skymax{i}'] = row.getValue("MAX")
            resultDict[f'skymin{i}'] = row.getValue("MIN")
            resultDict[f'skyave{i}'] = row.getValue("MEAN")
            clear_memory([row,rows])

        # Calculate zonal stats
        print(f'{PREFIX}Calculating mean anthropogenic sky luminance by zone...')
        outputTable = f"{gridsetp}skyzones.dbf"
        _ = arcpy.sa.ZonalStatisticsAsTable(
            zoneFile,"band",anthRaster,outputTable,"DATA","MEAN"
        )

        # Extract stats from output table
        zoneAve,zoneMax = [],[]
        rows = arcpy.SearchCursor(outputTable)
        for i in range(5):
            row = rows.next()
            zoneAve.append(row.getValue("MEAN"))
            zoneMax.append(row.getValue("COUNT"))
        clear_memory([row,rows])

        # Calculate zonal stats for each mosaic
        print(f'{PREFIX}Calculating all-sky anthropogenic luminous emittance...')
        mlux = nl_to_mlux(anthRaster)
        mlux70 = nl_to_mlux(anthRaster70)
        mlux80 = nl_to_mlux(anthRaster80)
        mosaicList = [mlux, mlux80, mlux70]
        for i,mosaic in enumerate(mosaicList):

            # Calculate zonal stats
            outputTable = f"{gridsetp}skyhemis{i}.dbf"
            _ = arcpy.sa.ZonalStatisticsAsTable(
                maskRaster,"VALUE",mosaic,outputTable,"DATA","ALL"
            )

            # Extract stats from output table
            rows = arcpy.SearchCursor(outputTable)
            row = rows.next()
            resultDict[f'hemis{i}'] = row.getValue("MEAN")
            resultDict[f'totalill{i}'] = row.getValue("SUM")
            clear_memory([row,rows])

        
        print(resultDict)
        print(zoneAve)
        print(zoneMax)
            
