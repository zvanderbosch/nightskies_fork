#-----------------------------------------------------------------------------#
#fullmosaic.py
#
#NPS Night Skies Program
#
#Last updated: 2025/07/24
#
#This script makes the whole sky mosaic from the full resolution images 
#according to the location of observed sky. The temporary files generated during
#the processing stage are stored in filepath.rasters/scratch_fullres folder. The 
#output raster and layer files are stored in file.griddata+dnight.
#
#
#Input: 
#   (1) ib###.tif, ib###.tfw
#           Full resolution TIFF files and associated TIFF World files (.tfw)
#           (filepath.calibdata/DATANIGHT/S_0#/tiff)
#   (2) pointerr_<DATASET>.txt
#           Pointing error data for image coordinates
#           (filepath.calibdata/DATANIGHT)
#   (3) extinction_fit_<FILTER>.xlsx
#           Best-fit extinction parameters
#           (filepath.calibdata/DATANIGHT)
#
#Output:
#   (1) skytopomags
#           Full resolution mosaic
#           (filepath.griddata/DATANIGHT/S_0#/fullres)
#   (2) skytopomags<DATASET>.lyrx 
#           Layer file for full-resolution mosaic
#           (filepath.griddata/DATANIGHT)
#
#History:
#	Dan Duriscoe -- Created as a module in firstbatchv4vb.py
#	Li-Wei Hung -- Cleaned and improved the code
#   Zach Vanderbosch -- Updated to Python 3.11 and ArcGIS Pro 3.5.2
#
#-----------------------------------------------------------------------------#

from glob import glob
from PIL import Image

import os
import stat
import arcpy
import shutil
import numpy as n
import pandas as pd

# Local Source
import filepath
import printcolors as pc

# Print status prefix
scriptName = 'fullmosaic.py'
PREFIX = f'{pc.GREEN}{scriptName:19s}{pc.END}: '

#-----------------------------------------------------------------------------#
if not os.path.exists(filepath.rasters+'scratch_fullres/'):
    os.makedirs(filepath.rasters+'scratch_fullres/')
    
# The geographic coordinate system WKT string
geogcs = (
    "GEOGCS["
        "'GCS_Sphere_EMEP',"
        "DATUM['D_Sphere_EMEP',"
        "SPHEROID['Sphere_EMEP',6370000.0,0.0]],"
        "PRIMEM['Greenwich',0.0],"
        "UNIT['Degree',0.0174532925199433]"
    "]"
)
          
#set arcpy environment variables part 1/2
arcpy.env.rasterStatistics = "NONE"
arcpy.env.overwriteOutput = True
arcpy.env.pyramid = "NONE"
arcpy.env.compression = "NONE"
arcpy.CheckOutExtension("3D")

# define source control points (37 points total, units = meters)
########################
#           *          #
#   *       *       *  #
#     *     *     *    #
#       *   *   *      #
#         * * *        #
# * * * * * * * * * * *#
#         * * *        #
#       *   *   *      #
#     *     *     *    #
#   *       *       *  #
#           *          #
########################
source_pnt = (
    "'0 0';"
    "'0 296039.8';'0 590759.1';'0 884157.9';'0 1176236';'0 1466994';"
    "'0 -296039.8';'0 -590759.1';'0 -884157.9';'0 -1176236';'0 -1466994';"
    "'-296039.8 0';'-590759.1 0';'-884157.9 0';'-1176236 0';'-1466994 0';"
    "'296039.8 0';'590759.1 0';'884157.9 0';'1176236 0';'1466994 0';"
    "'1241985 1241985';'-1241985 -1241985';'-1241985 1241985';'1241985 -1241985';"
    "'1445714 1445714';'-1445714 1445714';'-1445714 -1445714';'1445714 -1445714';"
    "'1037322 1037322';'-1037322 1037322';'-1037322 -1037322';'1037322 -1037322';"
    "'417730 417730';'-417730 417730';'-417730 -417730';'417730 -417730'"
)

# define target control points (37 points total)
target_pnt = (
    "'0 0';"
    "'0 296708';'0 593400';'0 890100';'0 1186800';'0 1483500';"
    "'0 -296700';'0 -593400';'0 -890100';'0 -1186800';'0 -1483500';"
    "'-296700 0';'-593400 0';'-890100 0';'-1186800 0';'-1483500 0';"
    "'296700 0';'593400 0';'890100 0';'1186800 0';'1483500 0';"
    "'1258791 1258791';'-1258791 -1258791';'-1258791 1258791';'1258791 -1258791';"
    "'1468590 1468590';'-1468590 1468590';'-1468590 -1468590';'1468590 -1468590';"
    "'1048993 1048993';'-1048993 1048993';'-1048993 -1048993';'1048993 -1048993';"
    "'419597 419597';'-419597 419597';'-419597 -419597';'419597 -419597'"
)
          
#------------------------------------------------------------------------------#
#-------------------            Define Functions            -------------------#
#------------------------------------------------------------------------------#

def clip_envelope(AZ, ALT, i):
    '''
    Function to generate rectangular clipping boundaries for an image.

    Parameters:
    -----------
    AZ: array
        Array of azimuth coordinates for each image center
    ALT: array
        Array of altitude coordinates for each image center
    i: int
        Image number

    Returns:
    --------
    boundary: str
        Coordinate boundaries for a given image
    '''
    if i < 15:
        bond = [AZ[i]-13,-6,AZ[i]+13,ALT[i]+12.9] 
    elif i < 30: 
        bond = [AZ[i]-13,ALT[i]-12.7,AZ[i]+13,ALT[i]+12.6] 
    elif i < 40:
        bond = [AZ[i]-18.6,ALT[i]-12.7,AZ[i]+18.6,ALT[i]+12.6] 
    elif i < 45:
        bond = [AZ[i]-39.6,ALT[i]-12.7,AZ[i]+39.6,ALT[i]+12.7] 
    boundary = ' '.join(str(i) for i in bond) #'xmin ymin xmax ymax'
    return boundary 
    
    
def tc(lon,lat):
    '''
    Returns the topocentric coordinate setting in WKT format

    Parameters:
    -----------
    lon: float
        Longitude coordinate
    lat: float
        Latitude coordinate

    Returns:
    --------
    topoCoord: str
        WKT formatted coordinate setting
    '''
    topoCoord = (
        "PROJCS["
            "'gnomonic',"
            f"{geogcs},"
            "PROJECTION['Gnomonic'],"
            "PARAMETER['False_Easting',0.0],"
            "PARAMETER['False_Northing',0.0],"
            f"PARAMETER['Longitude_Of_Center',{str(lon)}],"
            f"PARAMETER['Latitude_Of_Center',{str(lat)}],"
            "UNIT['Meter',1.0]"
        "]"
    )
    return topoCoord


def set_null_values(rasterFile):
    '''
    Function to set values within raster File to NoData

    Parameters:
    -----------
    rasterFile: string
        Path to raster File
    '''
    outSetNull = arcpy.sa.SetNull(
        rasterFile, 
        rasterFile, 
        "VALUE <= 0"
    )
    outSetNull.save(rasterFile)


def remove_readonly(func, path, excinfo):
    '''
    Error-catching function to handle removal of read-only folders

    Parameters:
    -----------
    func: python function
        Function to execute on path after chmod operation
    path: str
        Path to operate on
    excinfo: unknown
        Unused, but required by shutil
    '''
    os.chmod(path, stat.S_IWRITE)
    func(path)
  
  
def clear_dir(dir_path):
    '''
    Function to clear out all files and folders from
    the specified directory.

    Parameters:
    -----------
    dir_path: str
        Directory path to operate on
    '''
    for root, dirs, files in os.walk(dir_path, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.chmod(os.path.join(root, name), stat.S_IWRITE)
            os.rmdir(os.path.join(root, name))
    

#------------------------------------------------------------------------------#
#-------------------              Main Program              -------------------#
#------------------------------------------------------------------------------#

def mosaic(dnight, sets, filter):
    '''
    This module creates the mosaic of full-resolution images for each data set.

    Parameters:
    -----------
    dnight: str
        Name of data night to process
    sets: list
        List of data sets to process
    filter: str
        Name of photometric filter
    '''
    #set arcpy environment variables part 2/2
    arcpy.CheckOutExtension("Spatial")
    arcpy.env.workspace = filepath.rasters+'scratch_fullres/'
    arcpy.env.scratchWorkspace = filepath.rasters+'scratch_fullres'

    #filter paths
    F = {'V':'', 'B':'B/'}
    f = {'V':'', 'B':'b'}

    # Read in the best-fit enxtinction parameters
    extinctionFile = f"{filepath.calibdata}{dnight}/extinction_fit_{filter}.xlsx"
    extinctionData = pd.read_excel(extinctionFile)

    # Clear out and/or create domain shapefile directories
    for s in sets:
        domainsetp = f"{filepath.calibdata}{dnight}/S_0{s[0]}/{F[filter]}/domains/"
        if os.path.exists(domainsetp):
            shutil.rmtree(domainsetp, onerror=remove_readonly)
        os.makedirs(domainsetp)
    
    # Iterate over each data set
    for s in sets:

        # Convert set number to integer
        setnum = int(s[0])

        # Define file/folder paths
        calsetp = f"{filepath.calibdata}{dnight}/S_{setnum:02d}/{F[filter]}"
        gridsetp = f"{filepath.griddata}{dnight}/S_{setnum:02d}/{F[filter]}fullres/"
        scratchsetp = f"{filepath.rasters}scratch_fullres/"
        domainsetp = f"{calsetp}/domains/"

        # Remove and/or create gridsetp directory
        if os.path.exists(gridsetp):
            shutil.rmtree(gridsetp, onerror=remove_readonly)
        os.makedirs(gridsetp)

        # Clear out the scratch directory
        clear_dir(scratchsetp)
                
        # Read in the registered images coordinates
        file = f"{filepath.calibdata}{dnight}/pointerr_{setnum}.txt"
        Obs_AZ, Obs_ALT = n.loadtxt(file, usecols=(3,4)).T
        Obs_AZ[n.where(Obs_AZ>180)] -= 360
        Obs_AZ[35] %= 360
        imnum = len(Obs_AZ)

        # Get zeropoint, platescale, and exposure time for data set
        zeropoint = extinctionData['zeropoint_default'].iloc[setnum-1]
        platescale = extinctionData['avg_scale'].iloc[setnum-1]
        exptime = n.float64(extinctionData['exptime'].iloc[setnum-1])
        
        # Status Update
        print(f"{PREFIX}Generating fullres rasters for {filter}-Band Set {setnum}...")

        #loop through each file in the set
        for w in range(imnum+1):

            v = w+1
            if w == 45:
                w = 35
                Obs_AZ[w] -= 360

            # Copy TIFF file to scratch directory
            arcpy.management.CopyRaster(
                f'{calsetp}tiff/ib{w+1:03d}.tif', 
                f'ib{v:03d}.tif',
                "DEFAULTS",
                "","","","",
                "16_BIT_UNSIGNED"
            )
            
            # Re-define projection to topocentric coordinates
            arcpy.management.DefineProjection(
                f'ib{v:03d}.tif',
                tc(Obs_AZ[w],Obs_ALT[w])
            )
            
            # Warp image to remove barrel distortion image
            arcpy.management.Warp(
                f'ib{v:03d}.tif', 
                source_pnt, 
                target_pnt, 
                f'ibw{v:03d}.tif', 
                "POLYORDER3", 
                "BILINEAR"
            )
            set_null_values(f'ibw{v:03d}.tif')

            # Reproject into GCS
            arcpy.management.ProjectRaster(
                f'ibw{v:03d}.tif', 
                f'fwib{v:03d}.tif', 
                geogcs, 
                "BILINEAR", 
                "0.0261"
            )
            set_null_values(f'fwib{v:03d}.tif')

            # Create a raster border for clipping
            os.makedirs(f'{domainsetp}ib{v:03d}/')
            arcpy.ddd.RasterDomain(
                f'fwib{v:03d}.tif',
                f'{domainsetp}ib{v:03d}/ib{v:03d}_domain',
                'POLYGON'
            )
            arcpy.management.EliminatePolygonPart(
                f'{domainsetp}ib{v:03d}/ib{v:03d}_domain',
                f'{domainsetp}ib{v:03d}/ib{v:03d}_border',
                "PERCENT",
                "",
                "10",
                "CONTAINED_ONLY"
            )
                                       
            # Clip to image boundary
            # rectangle = clip_envelope(Obs_AZ, Obs_ALT, w)
            # arcpy.management.Clip("fwib%03d.tif"%v, rectangle, "fcib%03d"%v)
            clipFile = f'{domainsetp}ib{v:03d}/ib{v:03d}_border'
            arcpy.management.Clip(
                f"fwib{v:03d}.tif", 
                "", 
                f"fcib{v:03d}", 
                clipFile,
                "0",
                "ClippingGeometry",
                "NO_MAINTAIN_EXTENT"
            )

            # Status update
            if (v == w+1) & (v % 5 == 0):
                print(f'{PREFIX}{filter}-Band Set {setnum}, {v}/{imnum} rasters complete')
            
        # Mosaic raster list must start with an image with max pixel value > 256
        v=1; mstart=1
        while v < (len(Obs_AZ)+1):
            tiff = Image.open(filepath.rasters+'scratch_fullres/ib%03d.tif' %v)
            im = n.array(tiff)
            if n.max(im) > 255:
                mstart = v
                break
            v+=1
                        
        # Mosaic raster list
        R1 = ';'.join(['fcib%03d' %i for i in range(mstart,47)])
        R2 = ';'.join(['fcib%03d' %i for i in range(1,mstart)])
        R = R1+';'+R2

        # Mosaic to topocentric coordinate image; save in Griddata\
        print(f"{PREFIX}Mosaicking fullres rasters for {filter}-Band Set {setnum}...")
        arcpy.management.MosaicToNewRaster(
            R, gridsetp, 'skytopo', geogcs, 
            "32_BIT_FLOAT", "0.0261", "1", "BLEND", "FIRST"
        )

        # Crop lower latitude limit to -6.0 degrees
        arcpy.management.Clip(
            f'{gridsetp}skytopo', 
            "-180.0 -6.0 180.0 90.0", 
            f'{gridsetp}skytopoc'
        )
        
        # Convert to magnitudes per square arc second
        print(f"{PREFIX}Converting mosaic to mag/arcsec^2 for {filter}-Band Set {setnum}...")
        psa = 2.5*n.log10((platescale*60)**2) # platescale adjustment
        stm1 = arcpy.sa.Raster(gridsetp + os.sep + 'skytopoc')
        stm2 = 2.5 * arcpy.sa.Log10(stm1 / exptime)
        skytopomags = zeropoint + psa - stm2
        skytopomags.save(gridsetp + os.sep + 'skytopomags')

        # Save mags mosaic to disk
        print(f"{PREFIX}Creating fullres mosaic layer file for {filter}-Band Set {setnum}...")
        layerName = f"{dnight}_{setnum}_fullres{f[filter]}"
        layerfile = f"{filepath.griddata}{dnight}/skytopomags{f[filter]}{setnum}.lyrx"
        symbologyLayer = f"{filepath.rasters}magnitudes.lyrx"
        arcpy.management.MakeRasterLayer(gridsetp+'skytopomags', layerName)
        arcpy.management.ApplySymbologyFromLayer(layerName, symbologyLayer)
        arcpy.management.SaveToLayerFile(layerName, layerfile, "RELATIVE")

        # Remove intermediate raster directories & files
        del stm1, stm2, skytopomags
        shutil.rmtree(f'{gridsetp}skytopo', onerror=remove_readonly)
        shutil.rmtree(f'{gridsetp}skytopoc', onerror=remove_readonly)
        os.remove(f'{gridsetp}skytopo.aux.xml')
        os.remove(f'{gridsetp}skytopoc.aux.xml')

        # Clear scratch directory again
        clear_dir(scratchsetp)

        # Status update
        print(f"{PREFIX}{filter}-Band Set {setnum} fullres mosaic {pc.CYAN}COMPLETE{pc.END}")

    
if __name__ == "__main__":
    mosaic('FCNA160803', ['1st',],'V')
    pass
