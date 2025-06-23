#-----------------------------------------------------------------------------#
#starsvis.py
#
#NPS Night Skies Program
#
#Last updated: 2025/05/16
#
#This script computes the number/fraction of stars visible.
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

from astropy.time import Time
from astropy.io import fits

import os
import stat
import arcpy
import pandas as pd
import astropy.units as u
import astropy.coordinates as coord

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
scriptName = 'starsvis.py'
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


def get_zenith_coords(imageFile):
    '''
    Function to compute zenith RA and Dec coordinates
    at the time and location of a given image.
    '''

    # Load image header
    H = fits.getheader(imageFile)

    # Set observing time and site
    obsTime = Time(H['JD'], format='jd', scale='utc')
    obsSite = coord.EarthLocation.from_geodetic(
        lon = H['LONGITUD']*u.deg,
        lat = H['LATITUDE']*u.deg,
        height = H['ELEVATIO']*u.m
    )

    # Define Alt-Az coordinate object at zenith
    zenithTopoCoord = coord.SkyCoord(
        az=0.0*u.deg,
        alt=90.0*u.deg,
        obstime=obsTime,
        location=obsSite,
        frame='altaz'
    )

    # Transform to ICRS reference frame
    zenithICRSCoord = zenithTopoCoord.transform_to(coord.ICRS())
    raZenith = zenithICRSCoord.ra.deg
    decZenith = zenithICRSCoord.dec.deg

    # Wrap the RA coordinate at 180-deg
    if raZenith > 180:
        raZenith = -(raZenith - 360)
    else:
        raZenith = -raZenith

    return raZenith,decZenith


def get_coordinate_system(centralMeridian=180.0, latitudeOfOrigin=90.0):
    '''
    Function to define an ArcGIS coordinate system in WKT format
    '''

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
            f"PARAMETER['Central_Meridian',{centralMeridian:.8f}],"
            f"PARAMETER['Latitude_Of_Origin',{latitudeOfOrigin:.8f}],"
            "UNIT['Meter',1.0]"
        "]"
    )

    return coordinateSystem


def background_equation(varName):
    '''
    Function for generating background estimation equations 
    in string format with Python syntax for input to 
    arcpy.management.CalculateField
    '''

    stringEq = (
        f'( (0.00131744 * (!{varName}! ** 4))) - '
        f'(0.0925575 * (!{varName}! ** 3)) + '
        f'(2.38963265 * (!{varName}! ** 2)) - '
        f'(26.2432168 * !{varName}!) + '
        f'104.91069174'
    )

    return stringEq


#------------------------------------------------------------------------------#
#-------------------              Main Program              -------------------#
#------------------------------------------------------------------------------#

def calculate_stars_visible(dnight,sets,filter):
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

    # Load in the mask and airmass rasters
    maskRaster = arcpy.sa.Raster(f"{filepath.griddata}{dnight}/mask/maskd.tif")
    airmassRaster = arcpy.sa.Raster(f"{filepath.rasters}airmassf")

    # Convert mask raster to shape file
    maskShape = "mask.shp"
    arcpy.conversion.RasterToPolygon(
        maskRaster, maskShape, "SIMPLIFY", "VALUE"
    )

    # Define paths to additional shape files
    starShape = f"{filepath.rasters}shapefiles/SAOJ200079.shp"
    flatmaskShape = f"{filepath.rasters}shapefiles/allskyf.shp"

    # Load in bestfit extinction parameters
    extinctionFile = f"{filepath.calibdata}{dnight}/extinction_fit_{filter}.xlsx"
    extinctionData = pd.read_excel(extinctionFile)

    # Loop through each data set
    numstarsOutput = []
    for s in sets:

        # Set path for grid datasets
        setnum = int(s[0])
        calsetp = f"{filepath.calibdata}{dnight}/S_{setnum:02d}/{F[filter]}"
        gridsetp = f"{filepath.griddata}{dnight}/S_{setnum:02d}/{F[filter]}"

        # Initial status update for data set
        print(f'{PREFIX}Processing {dnight} Set-{setnum} {filter}-band...')

        # Get Zenith RA and Dec at dataset midpoint in time
        midpointImage = f"{calsetp}ib022.fit"
        raZenith, decZenith = get_zenith_coords(midpointImage)
        
        # Get extinction coefficient
        extCoeff = abs(extinctionData['extinction_fixedZ'].iloc[setnum-1])

        # Load in median sky brightness and natural sky rasters
        brightRasterMag = arcpy.sa.Raster(f"{gridsetp}median/skybrightmags")
        natRasterMag = arcpy.sa.Raster(f"{gridsetp}nat/natskymags")
        
        # Get coordinate systems
        coordSysLocal = get_coordinate_system()
        coordSysZenith = get_coordinate_system(
            centralMeridian=raZenith,
            latitudeOfOrigin=decZenith
        )

        # Project sky brightness rasters to fisheye equal area
        brightMagName = f"brightf{setnum}"
        arcpy.management.ProjectRaster(
            brightRasterMag, brightMagName, coordSysLocal, "BILINEAR", "5558.8"
        )
        brightRasterMagProjected = arcpy.sa.Raster(brightMagName)

        # Project star shapefile to local horizon
        starName = "j2000proj.shp"
        arcpy.management.Project(
            starShape, starName, coordSysZenith
        )
        arcpy.management.DefineProjection(
            starName, coordSysLocal
        )

        # Clip projected stars to flat and true horizons
        skyStarsFlatFile = "skystarsflat.shp"
        skyStarsFile = f"{gridsetp}nat/skystars.shp"
        arcpy.analysis.Clip(
            starName, flatmaskShape, skyStarsFlatFile
        )
        arcpy.analysis.Clip(
            starName, maskShape, skyStarsFile
        )
        numstarsf = int(arcpy.management.GetCount(skyStarsFlatFile).getOutput(0)) # Num stars to flat horizon
        numstarsm = int(arcpy.management.GetCount(skyStarsFile).getOutput(0))     # Num stars to observed horizon

        # Calculate extinction of stars
        arcpy.sa.ExtractMultiValuesToPoints( # Extract airmass values
            skyStarsFile, [[airmassRaster, "airmass"]], "NONE")
        arcpy.management.AddField(           # Add extmag field
            skyStarsFile, "extmag", "DOUBLE")
        arcpy.management.CalculateField(     # Calculate extincted star magnitudes
            skyStarsFile, "extmag", f"(!airmass! * {extCoeff:.8f}) + !vmag!", "PYTHON")

        # Extract background brightness values
        arcpy.sa.ExtractMultiValuesToPoints(
            skyStarsFile, [[brightRasterMagProjected, "background"]], "NONE")
        arcpy.sa.ExtractMultiValuesToPoints(
            skyStarsFile, [[natRasterMag, "backn"]], "NONE")

        # Calculate limiting magnitude at each star location
        arcpy.management.AddField(
            skyStarsFile, "lm", "DOUBLE")
        arcpy.management.CalculateField(
            skyStarsFile, "lm", background_equation('background'), "PYTHON")
        arcpy.management.AddField(
            skyStarsFile, "lmn", "DOUBLE")
        arcpy.management.CalculateField(
            skyStarsFile, "lmn", background_equation('backn'), "PYTHON")

        # Calculate visibility of each star
        arcpy.management.AddField(
            skyStarsFile, "em", "DOUBLE")
        arcpy.management.CalculateField(
            skyStarsFile, "em", '!lm! - !extmag!', "PYTHON")
        arcpy.management.AddField(
            skyStarsFile, "emn", "DOUBLE")
        arcpy.management.CalculateField(
            skyStarsFile, "emn", '!lmn! - !extmag!', "PYTHON")
        
        # Make a layer from the feature class
        skyStarsLayer = "skystars_lyr"
        arcpy.management.MakeFeatureLayer(
            skyStarsFile, skyStarsLayer)
        arcpy.management.SaveToLayerFile(
            skyStarsLayer, f"{gridsetp}nat/skystarslyr.lyrx", "ABSOLUTE")
        
        # Select stars with extincted mag < 7.5 only
        extStarsShape = f"{gridsetp}nat/extstars.shp"
        extStarsLayerName = "extstars_lyr"
        extStarsLayerFile = f"{gridsetp}nat/extstarslyr.lyrx"
        arcpy.management.SelectLayerByAttribute(
            skyStarsLayer, "NEW_SELECTION", "extmag < 7.9") # Is this a Bug?? Should be 7.5??
        arcpy.management.CopyFeatures(
            skyStarsLayer, extStarsShape)
        arcpy.management.MakeFeatureLayer(
            extStarsShape, extStarsLayerName)
        arcpy.management.SaveToLayerFile(
            extStarsLayerName, extStarsLayerFile, "ABSOLUTE")
        numstarse = int(arcpy.management.GetCount(extStarsShape).getOutput(0))

        # Select visible stars only
        visStarsShape = f"{gridsetp}nat/visstars.shp"
        visStarsLayerName = "visstars_lyr"
        visStarsLayerFile = f"{gridsetp}nat/visstarslyr.lyrx"
        arcpy.management.SelectLayerByAttribute(
            skyStarsLayer, "CLEAR_SELECTION", "")
        arcpy.management.SelectLayerByAttribute(
            skyStarsLayer, "NEW_SELECTION", "em > 0")
        arcpy.management.CopyFeatures(
            skyStarsLayer, visStarsShape)
        arcpy.management.MakeFeatureLayer(
            visStarsShape, visStarsLayerName)
        arcpy.management.SaveToLayerFile(
            visStarsLayerName, visStarsLayerFile, "ABSOLUTE")
        numstars = int(arcpy.GetCount_management(visStarsLayerFile).getOutput(0))

        # Get number of visible stars in natural sky
        natStarsShape = f"{gridsetp}nat/visstarsn.shp"
        natStarsLayerName = "visstarsn_lyr"
        natStarsLayerFile = f"{gridsetp}nat/visstarsnlyr.lyrx"
        arcpy.management.SelectLayerByAttribute(
            skyStarsLayer, "NEW_SELECTION", "emn > 0")
        arcpy.management.CopyFeatures(
            skyStarsLayer, natStarsShape)
        arcpy.management.MakeFeatureLayer(
            natStarsShape, natStarsLayerName)
        arcpy.management.SaveToLayerFile(
            natStarsLayerName, natStarsLayerFile, "ABSOLUTE")
        numstarsn = int(arcpy.GetCount_management(natStarsLayerFile).getOutput(0))

        # Apply symbologies to visible star layers
        symbologyVisStars = f"{filepath.rasters}skystarslp.lyrx"
        symbologyNatStars = f"{filepath.rasters}skystarslpn.lyrx"
        arcpy.ApplySymbologyFromLayer_management (visStarsLayerFile, symbologyVisStars)
        arcpy.ApplySymbologyFromLayer_management (natStarsLayerFile, symbologyNatStars)

        # Print out results
        print(f'{PREFIX}Number of stars to flat horizon            : {numstarsf}')
        print(f'{PREFIX}Number of stars to observed horizon        : {numstarsm}')
        print(f'{PREFIX}Number of stars visible without background : {numstarse}')
        print(f'{PREFIX}Number of stars visible in polluted sky    : {numstars}')
        print(f'{PREFIX}Number of stars visible in natural sky     : {numstarsn}')
        print(f'{PREFIX}Percent of stars visible in polluted sky   : {numstars/numstarsn*100:.1f} %')

        # Generate dataframe entry for given dataset
        numstarsEntry = pd.DataFrame(
            {
                'datanight': dnight,
                'dataset': setnum,
                'filter': filter,
                'Nstar_flat_horizon': numstarsf,
                'Nstar_obs_horizon': numstarsm,
                'Nstar_vis_noBkg': numstarse,
                'Nstar_vis_polluted': numstars,
                'Nstar_vis_natsky': numstarsn,
                'Nstar_vis_fraction': numstars / numstarsn
            },
            index = [setnum-1]
        )
        numstarsOutput.append(numstarsEntry)

    # Create final dataframe output
    numstarsOutput = pd.concat(numstarsOutput)

    return numstarsOutput