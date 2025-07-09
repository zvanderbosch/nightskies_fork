#-----------------------------------------------------------------------------#
#drawmaps.py
#
#NPS Night Skies Program
#
#Last updated: 2025/05/20
#
#This script generates the final output graphics including
#all-sky panaoramas and vertical illuminance figures.
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
#	Zach Vanderbosch -- Created script (translated from secondbatchv4.py
#                       and makepanoramas.vbs)
#
#-----------------------------------------------------------------------------#

from astropy.time import Time
from astropy.io import fits
from scipy.interpolate import make_interp_spline

import os
import sys
import stat
import arcpy
import numpy as n
import pandas as pd
import astropy.units as u
import matplotlib.pyplot as plt

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
scriptName = 'drawmaps.py'
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


def get_coordinate_system(centAz):
    '''
    Function to define an ArcGIS coordinate system in WKT format

    Parameters
    ----------
    centAz: int
        Central Azimuth of panorama, must be in range [0,360]
    '''

    # Check if centAz in allowed range
    if (centAz < 0) | (centAz > 360):
        print(f'{PREFIX}Central Azimuth {centAz} outside allowed range of 0 to 360')
        sys.exit(1)

    # Generate inputs to coordinate string
    if centAz <= 180:
        sysName = f"HA {centAz:03d}"
        centralMeridian = float(centAz)
    elif (centAz > 180) & (centAz < 360):
        sysName = f"HA {centAz:03d}"
        centralMeridian = float(centAz - 360)
    elif centAz == 360:
        sysName = f"HA {centAz:03d}"
        centralMeridian = 0.0

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
            f"'{sysName}',"
            f"{geogcs},"
            "PROJECTION['Hammer_Aitoff'],"
            "PARAMETER['False_Easting',0.0],"
            "PARAMETER['False_Northing',0.0],"
            f"PARAMETER['Central_Meridian',{centralMeridian:.1f}],"
            "UNIT['Meter',1.0]"
        "]"
    )

    return coordinateSystem


def get_site_info(imageFile):
    '''
    Function that computes the Local Mean Time (LMT)
    and returns the site location name and observer
    names along with a formatted string containing 
    date and time info.
    '''

    # Load image header
    H = fits.getheader(imageFile)

    # Get site and observer names
    siteName = H['LOCATION']
    observers = H['OBSERVER']

    # Set observing time
    t = Time(
        H['DATE-OBS'], 
        format='isot', 
        scale='utc'#,
    )

    # Convert UTC to LMT using site longitude
    dayShift = 0
    hourfracUTC = t.ymdhms.hour + t.ymdhms.minute/60 + t.ymdhms.second/3600
    hourfracLMT = hourfracUTC + H['LONGITUD']/15
    if hourfracLMT < 0:
        hourfracLMT += 24
        dayShift -= 1

    # Approximate LMT midpoint of dataset
    midpointLMT = hourfracLMT + 0.16
    if midpointLMT > 24:
        midpointLMT -= 24
        dayShift += 1

    # Apply day shift
    t = t + dayShift*u.day

    # Generate the full date-time string
    dt = t.datetime
    day = dt.strftime("%d").lstrip('0')
    month = dt.strftime("%B")
    year = dt.strftime("%Y")
    dateString = f"{month} {day}, {year}  {midpointLMT:.1f} hours LMT"

    return siteName, observers, dateString


def make_vertillum_figure(dataNight, setNumber):
    '''
    Function to generate figure comparing
    vertical illuminance from all light sources
    and from anthropogenic light only.
    '''

    # Status update
    print(f'{PREFIX}Generating final vertical illummination figure...')

    # Load in vertical illuminance Excel sheet (vert.xlsx)
    vertFile = f"{filepath.calibdata}{dataNight}/vert.xlsx"
    vertDataAll = pd.read_excel(
        vertFile, 
        sheet_name="Allsky_All_Sources", 
        skiprows=10,
        header=None,
        names=['Azimuth','VertIllum'],
        usecols=(0,setNumber)
    )
    vertDataAnth = pd.read_excel(
        vertFile, 
        sheet_name="Allsky_Artificial", 
        skiprows=10,
        header=None,
        names=['Azimuth','VertIllum'],
        usecols=(0,setNumber)
    )

    # Copy values at 0 degrees to 360 degrees
    vertDataAll.loc[len(vertDataAll)] = [
        360., vertDataAll['VertIllum'].iloc[0]
    ]
    vertDataAnth.loc[len(vertDataAnth)] = [
        360., vertDataAnth['VertIllum'].iloc[0]
    ]

    # Get smoothed spline interpolations of data
    xsmooth = n.linspace(0.,360.,1000)
    splineInterpAll = make_interp_spline(
        vertDataAll['Azimuth'].values, 
        vertDataAll['VertIllum'].values
    )
    splineInterpAnth = make_interp_spline(
        vertDataAnth['Azimuth'].values, 
        vertDataAnth['VertIllum'].values
    )
    vertInterpAll = splineInterpAll(xsmooth)
    vertInterpAnth = splineInterpAnth(xsmooth)

    # Create figures to plot horizontal and vertical illuminances
    fig = plt.figure('vert',figsize=(10,6))
    ax = fig.add_subplot(111)

    # Plot smoothed vertical illuminance
    ax.plot(xsmooth, vertInterpAll , ls='-', lw=3, c='b', label='All Sources to 96 ZA')
    ax.plot(xsmooth, vertInterpAnth, ls='-', lw=3, c='r', label='Artificial Skyglow Only')

    # Add legend
    ax.legend(
        loc='upper right', 
        handlelength=2.0,
        fontsize=11, 
        ncol=1,
        facecolor='gainsboro',
        framealpha=1.0
    )

    # Adjust plot appearances
    ax.set_xlim(0, 360)
    ax.set_ylim(0, 1.15*ax.get_ylim()[1])
    ax.set_xticks(n.arange(0,361,30))
    ax.set_axisbelow(True)
    ax.grid(ls=':', c='silver', lw=0.75)
    ax.set_xlabel('Azimuth (deg)',fontsize=14)
    ax.set_ylabel('Vertical Illumination (mLux)',fontsize=14)
    ax.set_title(f"{dataNight} Data Set {setNumber}", fontsize=16)

    # Save figure to file
    fileName = f"{filepath.graphics}{dataNight}_vert_{setNumber}.png"
    plt.savefig(fileName, dpi=300, bbox_inches='tight')
    plt.close('vert')


def make_fullres_panorama(mapObj, layoutObj, subtitleTextElement, gridPath, dataNight, setNumber, centAz):
    '''
    Function to generate full resolution panorama graphic
    '''

    # Add full-resolution mosaic to template
    subtitleTextElement.text = "Full Resolution Mosaic"
    mosaicLayer = arcpy.mp.LayerFile(
        f"{gridPath}skytopomags{setNumber}.lyrx"
    )
    mapObj.addLayer(mosaicLayer, "BOTTOM")

    # Export template to JPEG
    print(f'{PREFIX}Generating JPEG for Full Resolution Mosaic...')
    jpegFile = f"{filepath.graphics}{dataNight}_fullres_{setNumber}_HA{centAz:03d}.jpg"
    layoutObj.exportToJPEG(
        jpegFile, resolution=200
    )

    # Remove the added layer
    for lyr in mapObj.listLayers():
        if "fullres" in lyr.name:
            mapObj.removeLayer(lyr)
    

def make_skybright_panorama(mapObj, layoutObj, subtitleTextElement, gridPath, dataNight, setNumber, centAz):
    '''
    Function to generate sky brightness panorama graphic
    '''

    # Add full-resolution mosaic to template
    subtitleTextElement.text = "Background Sky Brightness"
    mosaicLayer = arcpy.mp.LayerFile(
        f"{gridPath}skybrightmags{setNumber}.lyrx"
    )
    mapObj.addLayer(mosaicLayer, "BOTTOM")

    # Export template to JPEG
    print(f'{PREFIX}Generating JPEG for Sky Brightness Mosaic...')
    jpegFile = f"{filepath.graphics}{dataNight}_skybright_{setNumber}_HA{centAz:03d}.jpg"
    layoutObj.exportToJPEG(
        jpegFile, resolution=200
    )

    # Remove the added layer
    for lyr in mapObj.listLayers():
        if "median" in lyr.name:
            mapObj.removeLayer(lyr)


def make_anthropogenic_panorama(mapObj, layoutObj, subtitleTextElement, gridPath, dataNight, setNumber, centAz):
    '''
    Function to generate anthropogenic light panorama graphic
    '''

    # Add full-resolution mosaic to template
    subtitleTextElement.text = "Estimated Artificial Sky Glow"
    mosaicLayer = arcpy.mp.LayerFile(
        f"{gridPath}anthlightmags{setNumber}.lyrx"
    )
    mapObj.addLayer(mosaicLayer, "BOTTOM")

    # Export template to JPEG
    print(f'{PREFIX}Generating JPEG for Anthropogenic Mosaic...')
    jpegFile = f"{filepath.graphics}{dataNight}_artificial_{setNumber}_HA{centAz:03d}.jpg"
    layoutObj.exportToJPEG(
        jpegFile, resolution=200
    )

    # Remove the added layer
    for lyr in mapObj.listLayers():
        if "skyglow" in lyr.name:
            mapObj.removeLayer(lyr)


def make_naturalsky_panorama(mapObj, layoutObj, subtitleTextElement, gridPath, dataNight, setNumber, centAz):
    '''
    Function to generate natural sky model panorama graphic
    '''

    # Add full-resolution mosaic to template
    subtitleTextElement.text = "Modeled Natural Sky Background Brightness"
    mosaicLayer = arcpy.mp.LayerFile(
        f"{gridPath}natskymags{setNumber}.lyrx"
    )
    mapObj.addLayer(mosaicLayer, "BOTTOM")

    # Export template to JPEG
    print(f'{PREFIX}Generating JPEG for Natural Sky Mosaic...')
    jpegFile = f"{filepath.graphics}{dataNight}_natsky_{setNumber}_HA{centAz:03d}.jpg"
    layoutObj.exportToJPEG(
        jpegFile, resolution=200
    )

    # Remove the added layer
    for lyr in mapObj.listLayers():
        if "natsky" in lyr.name:
            mapObj.removeLayer(lyr)



#-----------------------------------------------------------------------------#
# Main Program

def generate_graphics(dnight,sets,processorName,centralAzimuth,locationName):
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

    # Define coordinate system based on central Azimuth
    centralAzimuth = int(n.floor(centralAzimuth))
    coordSys = get_coordinate_system(centralAzimuth)

    # Load in panorama template ArcGIS project
    templateMap = f"{filepath.maps}templatemap/templatemap.aprx"
    p = arcpy.mp.ArcGISProject(templateMap)

    # Set new spatial reference for the map
    mxd = p.listMaps("Layers")[0]
    newSR = arcpy.SpatialReference(text=coordSys)
    mxd.spatialReference = newSR

    # Set extent of the panorama
    layout = p.listLayouts()[0]
    mapFrame = layout.listElements("mapframe_element")[0]
    extent = mapFrame.camera.getExtent()
    extent.XMin = -18055000
    extent.XMax = 18055000
    extent.YMin = 30000
    extent.YMax = 8050000
    mapFrame.camera.setExtent(extent)

    # Grab individual map elements
    mapTitle = layout.listElements("text_element","Title")[0]
    mapSubtitle = layout.listElements("text_element","Subtitle")[0]
    mapCollBy = layout.listElements("text_element","collected by")[0]
    mapProcBy = layout.listElements("text_element","processed by")[0]


    # Loop through each data set
    for s in sets:

        # Set paths for calibrated images and grid datasets
        setnum = int(s[0])
        calsetp = f"{filepath.calibdata}{dnight}/S_{setnum:02d}/{F['V']}"
        gridsetp = f"{filepath.griddata}{dnight}/"

        # Initial status update for data set
        print(f'{PREFIX}Generating graphics for {dnight} Set-{setnum}...')

        # Get Zenith RA and Dec at dataset midpoint in time
        firstImage = f"{calsetp}zenith1.fit"
        site, observers, date = get_site_info(firstImage)

        # Generate map title
        titleText = f"{locationName.replace('_',' ')}  {site}  {date}"

        # Set text elements
        mapTitle.text = titleText
        mapCollBy.text = observers
        mapProcBy.text = processorName.replace("_"," ")

        # Make vertical illumination figure
        make_vertillum_figure(dnight, setnum)

        # Make full-resolution panorama
        make_fullres_panorama(mxd, layout, mapSubtitle, gridsetp, dnight, setnum, centralAzimuth)

        # Make sky-brightness panorama
        make_skybright_panorama(mxd, layout, mapSubtitle, gridsetp, dnight, setnum, centralAzimuth)

        # Make anthropogenic light panorama
        make_anthropogenic_panorama(mxd, layout, mapSubtitle, gridsetp, dnight, setnum, centralAzimuth)

        # Make natural sky model panorama
        make_naturalsky_panorama(mxd, layout, mapSubtitle, gridsetp, dnight, setnum, centralAzimuth)
