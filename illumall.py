#-----------------------------------------------------------------------------#
#skyglow.py
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

    # Set mosaic processing order
    mosaicOrder = ['allsky','za80','za70']

    # Iterate over each mosaic
    for i,key in enumerate(mosaicOrder):
        
        # Choose the mosaic to process
        mosaic = mosaicDict[key]

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
    Function for calculating luminous emittance metrics
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


def calc_horizontal_illuminance(mosaicDict, zoneRaster, gridPath, results):
    '''
    Function for calculating horizontal illuminance metrics
    '''

    # Load in horizontal illuminance factor grid
    horizRasterFile = f"{filepath.rasters}horizillf"
    horizRaster = arcpy.sa.Raster(horizRasterFile)

    # Convert from nL to mlux units and multiply by factor grid
    mlux = nl_to_mlux(mosaicDict['allsky'])
    mlux80 = nl_to_mlux(mosaicDict['za80'])
    mlux70 = nl_to_mlux(mosaicDict['za70'])
    mosaicListHoriz = [
        mlux * horizRaster, 
        mlux80 * horizRaster, 
        mlux70 * horizRaster
    ]

    # Iterate over each mosaic
    for i,mosaic in enumerate(mosaicListHoriz):

        # Calculate zonal stats
        outputTable = f"{gridPath}skyhemis{i}.dbf"
        _ = arcpy.sa.ZonalStatisticsAsTable(
            zoneRaster,"VALUE",mosaic,outputTable,"DATA","SUM"
        )

        # Extract stats from output table
        rows = arcpy.SearchCursor(outputTable)
        row = rows.next()
        results[f'horizs{i}'] = row.getValue("SUM")
        clear_memory([row,rows])

    return results


def calc_vertical_illuminance(mosaicDict, zoneRaster, gridPath, results):
    '''
    Function for calculating vertical illuminance metrics
    '''

    # Convert skyglow from nL to mlux units
    mlux = nl_to_mlux(mosaicDict['allsky'])
    mlux80 = nl_to_mlux(mosaicDict['za80'])
    mlux70 = nl_to_mlux(mosaicDict['za70'])

    # Define array of rotation angles (0 -> 350 in 10-degree increments)
    rotationAngles = n.arange(0,360,10,dtype=int)

    # Iterate over rotation angles in 10-degree increments
    for angle in rotationAngles:

        # Status update
        print(f'{PREFIX}Calculating vertical illuminance at angle {angle}...')

        # Load in corresponding raster file
        vertRasterFile = f"{filepath.rasters}vertgrids/vertf{-angle}"
        vertRaster = arcpy.sa.Raster(vertRasterFile)

        # Calculate vertical illuminance in mlux for each zone
        mosaicListVert = [
            mlux * vertRaster,
            mlux80 * vertRaster,
            mlux70 * vertRaster
        ]

        # Iterate over each mosaic
        for i,mosaic in enumerate(mosaicListVert):

            # Calculate zonal stats
            outputTable = f"{gridPath}skyvert{i}.dbf"
            _ = arcpy.sa.ZonalStatisticsAsTable(
                zoneRaster,"VALUE",mosaic,outputTable,"DATA","SUM"
            )

            # Extract stats from output table
            rows = arcpy.SearchCursor(outputTable)
            row = rows.next()
            results[f'vert-{angle:03d}-{i}'] = row.getValue("SUM")
            clear_memory([row,rows])

    return results


def calc_sqi_histograms(zoneRaster, gridPath, clipFile80, clipFile70):
    '''
    Function to calculate Sky Quality Index (SQI) histograms
    '''

    # Clip skyglow magnitudes mosaic to 80 and 70 Zenith Angle limits
    print(f'{PREFIX}Clipping skyglow (mag) mosaic to 80 and 70 Zenith Angle limits...')
    inputRaster = f"{gridPath}skyglow/anthlightmags"
    arcpy.management.Clip(
        inputRaster,"#", "anth80mags", clipFile80,"#","ClippingGeometry"
    )
    arcpy.management.Clip(
        inputRaster,"#", "anth70mags", clipFile70,"#","ClippingGeometry"
    )

    # Create raster layers from anthropogenic mags mosaic
    print(f'{PREFIX}Generating skyglow (mag) raster layers...')
    layerName = "anthmaskedlyr"
    layerName80 = "anth80lyr"
    layerName70 = "anth70lyr"
    symbologyLayer = f"{filepath.rasters}magnitudes1.lyrx"
    arcpy.management.MakeRasterLayer(
        inputRaster,layerName,"#","#","#"
    )
    arcpy.management.MakeRasterLayer(
        "anth80mags",layerName80,"#","#","#"
    )
    arcpy.management.MakeRasterLayer(
        "anth70mags",layerName70,"#","#","#"
    )
    arcpy.management.ApplySymbologyFromLayer(
        layerName, symbologyLayer
    )
    arcpy.management.ApplySymbologyFromLayer(
        layerName80, symbologyLayer
    )
    arcpy.management.ApplySymbologyFromLayer(
        layerName70, symbologyLayer
    )

    # Iterate over each layer
    print(f'{PREFIX}Calculating SQI zonal histograms...')
    layerList = [layerName, layerName80, layerName70]
    for layer in layerList:

        # Determine output table name
        if "80" in layer:
            outputTable = f"{gridPath}sqitbl80.dbf"
        elif "70" in layer:
            outputTable = f"{gridPath}sqitbl80.dbf"
        else:
            outputTable = f"{gridPath}sqitbl.dbf"

        # Generate zonal histogram
        arcpy.sa.ZonalHistogram(
            zoneRaster, "VALUE", layer, outputTable, ""
        )

        # Delete the layer
        arcpy.management.Delete(layer,"")


#-----------------------------------------------------------------------------#
# Main Program

def calculate_illuminance(dnight,sets,filter):
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
        print(f'{PREFIX}Clipping skyglow (nL) mosaic to 80 and 70 Zenith Angle limits...')
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
        resultDict = calc_sky_luminance(mosaicDict, maskRaster, gridsetp, resultDict)


        # Calculate mean sky luminance by zone
        print(f'{PREFIX}Calculating mean anthropogenic sky luminance by zone...')
        resultDict = calc_zonal_sky_luminance(anthRaster, zoneFile, gridsetp, resultDict)


        # Calculate zonal stats for each mosaic
        print(f'{PREFIX}Calculating all-sky anthropogenic luminous emittance...')
        resultDict = calc_luminouos_emittance(mosaicDict, maskRaster, gridsetp, resultDict)


        # Calculate horizontal illuminance
        print(f'{PREFIX}Calculating anthropogenic horizontal illuminance...')
        resultDict = calc_horizontal_illuminance(mosaicDict, maskRaster, gridsetp, resultDict)


        # Calculate horizontal illuminance
        print(f'{PREFIX}Calculating anthropogenic vertical illuminance...')
        resultDict = calc_vertical_illuminance(mosaicDict, maskRaster, gridsetp, resultDict)


        # Calculate zonal histograms
        print(f'{PREFIX}Calculating anthropogenic SQI zonal histograms...')
        calc_sqi_histograms(maskRaster, gridsetp, za80File, za70File)
            
