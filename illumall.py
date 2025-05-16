#-----------------------------------------------------------------------------#
#illumall.py
#
#NPS Night Skies Program
#
#Last updated: 2025/05/15
#
#This script computes sky luminance and illuminance
#statistics for all light sources using the skybright
#median filtered mosaic.
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

import os
import stat
import arcpy
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
scriptName = 'illumall.py'
PREFIX = f'{pc.GREEN}{scriptName:19s}{pc.END}: '

#------------------------------------------------------------------------------#
#-----------------  Define Fisheye Equal Area Coord System  -------------------#
#------------------------------------------------------------------------------#

geogcs = (
    "GEOGCS["
        "'GCS_Sphere_EMEP',"
        "DATUM['D_Sphere_EMEP',"
        "SPHEROID['Sphere_EMEP',6370000.0,0.0]],"
        "PRIMEM['Greenwich',0.0],"
        "UNIT['Degree',0.0174532925199433]"
    "]"
)
coordinateSystem = (
    "PROJCS["
        "'fisheye equal area',"
        f"{geogcs},"
        "PROJECTION['Lambert_Azimuthal_Equal_Area'],"
        "PARAMETER['False_Easting',0.0],"
        "PARAMETER['False_Northing',0.0],"
        "PARAMETER['Central_Meridian',180.0],"
        "PARAMETER['Latitude_Of_Origin',90.0],"
        "UNIT['Meter',1.0]"
    "]"
)

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


def nl_to_mlux(raster):
    '''
    Unit conversion from nano-Lamberts to milli-Lux
    '''
    return (0.000000761544 * raster) / 314.159


def calc_sky_area_luminance(mosaic, zoneRaster, gridPath, results):
    '''
    Function for calculating observed sky luminance is full area
    down to Zenith Angle 96, including below horizon values.
    '''

    # Calculate zonal stats
    outputTable = f"{gridPath}skyhemisall.dbf"
    _ = arcpy.sa.ZonalStatisticsAsTable(
        zoneRaster,"VALUE",mosaic,outputTable,"DATA","ALL"
    )

    # Extract stats from output table
    rows = arcpy.SearchCursor(outputTable)
    row = rows.next()
    results[f'skymax0'] = row.getValue("MAX")
    results[f'skymin0'] = row.getValue("MIN")
    results[f'skyave0'] = row.getValue("MEAN")
    clear_memory([row,rows])

    return results


def calc_sky_mask_luminance(mosaic, zoneRaster, gridPath, results):
    '''
    Function for calculating observed all-sky luminance,
    excluding below horizon values.
    '''

    # Calculate zonal stats
    outputTable = f"{gridPath}skyhemismasked.dbf"
    _ = arcpy.sa.ZonalStatisticsAsTable(
        zoneRaster,"VALUE",mosaic,outputTable,"DATA","ALL"
    )

    # Extract stats from output table
    rows = arcpy.SearchCursor(outputTable)
    row = rows.next()
    results[f'skymax1'] = row.getValue("MAX")
    results[f'skymin1'] = row.getValue("MIN")
    results[f'skyave1'] = row.getValue("MEAN")
    clear_memory([row,rows])

    return results


def calc_area_luminouos_emittance(mosaic, zoneRaster, gridPath, results):
    '''
    Function for calculating observed luminous emittance in full
    area down to Zenith Angle 96, including below horizon values.
    '''

    # Convert from nL to mlux units
    mosaicMlux = nl_to_mlux(mosaic)

    # Calculate zonal stats
    outputTable = f"{gridPath}skyhemisall1.dbf"
    _ = arcpy.sa.ZonalStatisticsAsTable(
        zoneRaster,"VALUE",mosaicMlux,outputTable,"DATA","ALL"
    )

    # Extract stats from output table
    rows = arcpy.SearchCursor(outputTable)
    row = rows.next()
    results[f'hemis0'] = row.getValue("MEAN")
    results[f'totalill0'] = row.getValue("SUM")
    clear_memory([row,rows])

    return results


def calc_mask_luminouos_emittance(mosaic, zoneRaster, gridPath, results):
    '''
    Function for calculating observed luminous emittance,
    excluding below-horizon values
    '''

    # Convert from nL to mlux units
    mosaicMlux = nl_to_mlux(mosaic)

    # Calculate zonal stats
    outputTable = f"{gridPath}skyhemisall2.dbf"
    _ = arcpy.sa.ZonalStatisticsAsTable(
        zoneRaster,"VALUE",mosaicMlux,outputTable,"DATA","ALL"
    )

    # Extract stats from output table
    rows = arcpy.SearchCursor(outputTable)
    row = rows.next()
    results[f'hemis1'] = row.getValue("MEAN")
    results[f'totalill1'] = row.getValue("SUM")
    clear_memory([row,rows])

    return results


def calc_scalar_illuminance(mosaic, zoneRaster, gridPath, results):
    '''
    Function for calculating observed scalar luminance, using
    all values down to Zenith Angle 90 (full hemisphere).
    '''

    # Convert from nL to mlux units
    mosaicMlux = nl_to_mlux(mosaic)

    # Calculate zonal stats
    outputTable = f"{gridPath}skyscalar.dbf"
    _ = arcpy.sa.ZonalStatisticsAsTable(
        zoneRaster,"VALUE",mosaicMlux,outputTable,"DATA","SUM"
    )

    # Extract stats from output table
    rows = arcpy.SearchCursor(outputTable)
    row = rows.next()
    results['skyscalar'] = row.getValue("SUM")
    clear_memory([row,rows])

    return results


def calc_horizontal_illuminance(mosaic, zoneRaster, gridPath, results):
    '''
    Function for calculating horizontal illuminance metrics
    '''

    # Load in horizontal illuminance factor grid
    horizRasterFile = f"{filepath.rasters}horizillf"
    horizRaster = arcpy.sa.Raster(horizRasterFile)

    # Convert from nL to mlux units and multiply by factor grid
    mosaicMlux = nl_to_mlux(mosaic)
    mosaicHoriz = mosaicMlux * horizRaster

    # Calculate zonal stats
    outputTable = f"{gridPath}skyhoriz.dbf"
    _ = arcpy.sa.ZonalStatisticsAsTable(
        zoneRaster,"VALUE",mosaicHoriz,outputTable,"DATA","SUM"
    )

    # Extract stats from output table
    rows = arcpy.SearchCursor(outputTable)
    row = rows.next()
    results['horizs'] = row.getValue("SUM")
    clear_memory([row,rows])

    return results


def calc_vertical_illuminance(mosaic, zoneRaster, gridPath, results):
    '''
    Function for calculating vertical illuminance metrics
    '''

    # Convert from nL to mlux units
    mosaicMlux = nl_to_mlux(mosaic)

    # Define array of rotation angles (0 -> 355 in 5-degree increments)
    rotationAngles = n.arange(0,360,5,dtype=int)

    # Iterate over rotation angles in 10-degree increments
    for angle in rotationAngles:

        # Status update
        print(f'{PREFIX}Calculating vertical illuminance at angle {angle}...')

        # Load in corresponding raster file
        vertRasterFile = f"{filepath.rasters}vertgrids2/vertf{-angle}"
        vertRaster = arcpy.sa.Raster(vertRasterFile)

        # Calculate vertical illuminance in mlux
        mosaicVert = mosaicMlux * vertRaster

        # Calculate zonal stats
        outputTable = f"{gridPath}skyvert0all.dbf"
        _ = arcpy.sa.ZonalStatisticsAsTable(
            zoneRaster,"VALUE",mosaicVert,outputTable,"DATA","SUM"
        )

        # Extract stats from output table
        rows = arcpy.SearchCursor(outputTable)
        row = rows.next()
        results[f'vert-{angle:03d}'] = row.getValue("SUM")
        clear_memory([row,rows])

    return results



#------------------------------------------------------------------------------#
#-------------------              Main Program              -------------------#
#------------------------------------------------------------------------------#

def calculate_statistics(dnight,sets,filter):
    '''
    Main program for computing sky luminance and illuminance
    statistics from only anthropogenic sources. 
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

    # Load in mask rasters for zonal stats
    maskRaster = arcpy.sa.Raster(f"{filepath.griddata}{dnight}/mask/maskd.tif")
    areaRaster = arcpy.sa.Raster(f"{filepath.rasters}arearasterf")
    hemiRaster = arcpy.sa.Raster(f"{filepath.rasters}hemirasterf")

    # Loop through each data set
    for s in sets:

        # Set path for grid datasets
        setnum = int(s[0])
        gridsetp = f"{filepath.griddata}{dnight}/S_{setnum:02d}/{F[filter]}"

        # Initial status update for data set
        print(f'{PREFIX}Processing {dnight} Set-{setnum} {filter}-band...')

        # Load in median sky brightness mosaic in nano-Lamberts [nL]
        inRaster = arcpy.sa.Raster(f"{gridsetp}median/skybrightnl")

        # Project raster into fisheye coordinate system
        brightRasterFile = "skynlf.tif"
        arcpy.ProjectRaster_management (
            inRaster, brightRasterFile, coordinateSystem, "BILINEAR","5558.8"
        )
        brightRaster = arcpy.sa.Raster(brightRasterFile)

        # Initialize result dict
        resultDict = {}

        # Calculate sky luminance metrics
        print(f'{PREFIX}Calculating sky luminance for all sources...')
        resultDict = calc_sky_area_luminance(brightRaster, areaRaster, gridsetp, resultDict)
        resultDict = calc_sky_mask_luminance(brightRaster, maskRaster, gridsetp, resultDict)


        # Calculate luminous emittance
        print(f'{PREFIX}Calculating luminous emittance for all sources...')
        resultDict = calc_area_luminouos_emittance(brightRaster, areaRaster, gridsetp, resultDict)
        resultDict = calc_mask_luminouos_emittance(brightRaster, maskRaster, gridsetp, resultDict)


        # Calculate scalar illuminance
        print(f'{PREFIX}Calculating scalar illuminance for all sources...')
        resultDict = calc_scalar_illuminance(brightRaster, hemiRaster, gridsetp, resultDict)


        # Calculate horizontal illuminance
        print(f'{PREFIX}Calculating horizontal illuminance for all sources...')
        resultDict = calc_horizontal_illuminance(brightRaster, hemiRaster, gridsetp, resultDict)


        # Calculate vertical illuminance
        print(f'{PREFIX}Calculating anthropogenic vertical illuminance...')
        resultDict = calc_vertical_illuminance(brightRaster, areaRaster, gridsetp, resultDict)
