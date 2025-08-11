#-----------------------------------------------------------------------------#
#medianmosaic.py
#
#NPS Night Skies Program
#
#Last updated: 2025/07/23
#
#This script makes the whole sky mosaic from the median filtered images 
#according to the location of observed sky. The temporary files generated during
#the processing stage are stored in filepath.rasters/scratch_median folder. The 
#output raster and layer files are stored in file.griddata+dnight.
#
#
#Input: 
#   (1) median_ib###.tif, median_ib###.tfw
#           Median-filtered TIFF files and associated TIFF World files (.tfw)
#           (filepath.calibdata/DATANIGHT/S_0#/tiff)
#   (2) pointerr_<DATASET>.txt
#           Pointing error data for image coordinates
#           (filepath.calibdata/DATANIGHT)
#   (3) extinction_fit_<FILTER>.xlsx
#           Best-fit extinction parameters
#           (filepath.calibdata/DATANIGHT)
#
#Output:
#   (1) skybright
#           Median-filtered sky brightness mosaic in ADU units after resampling
#           (filepath.griddata/DATANIGHT/S_0#/median)
#   (2) skybrightmags
#           Median-filtered sky brightness mosaic (mag/arcsec^2)
#           (filepath.griddata/DATANIGHT/S_0#/median)
#   (3) skybrightmags<DATASET>.lyrx 
#           Layer file for median-filtered sky brightness mosaic (mag/arcsec^2)
#           (filepath.griddata/DATANIGHT)
#   (4) mask.tif
#           Used to make horizon mask with PhotoShop in later process
#           (filepath.griddata/DATANIGHT/mask)
#
#
#History:
#	Dan Duriscoe -- Created as a module in firstbatchv4vb.py
#	Li-Wei Hung -- Cleaned and improved the code
#   Zach Vanderbosch -- Updated to Python 3.11 and ArcGIS Pro 3.5.2
#
#-----------------------------------------------------------------------------#
from glob import glob
from astropy.io import fits
from PIL import Image
from skimage.transform import downscale_local_mean

import os
import stat
import time
import arcpy
import shutil
import numpy as n
import pandas as pd

# Local Source
import filepath
import printcolors as pc

# Print status prefix
scriptName = 'medianmosaic.py'
PREFIX = f'{pc.GREEN}{scriptName:19s}{pc.END}: '

#-----------------------------------------------------------------------------#
if not os.path.exists(filepath.rasters+'scratch_median/'):
    os.makedirs(filepath.rasters+'scratch_median/')
    
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
arcpy.env.rasterStatistics = "STATISTICS 2 2 (-999)"
arcpy.env.overwriteOutput = True
arcpy.env.pyramid = "NONE"
arcpy.env.compression = "NONE"

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


def clear_scratch(scratch_dir):
    '''
    Function to clear out all files and folders from
    the scratch directory.

    Parameters:
    -----------
    scratch_dir: str
        Directory path to operate on
    '''
    for root, dirs, files in os.walk(scratch_dir, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.chmod(os.path.join(root, name), stat.S_IWRITE)
            os.rmdir(os.path.join(root, name))


def check_file(filename):
    ''' 
    Check that a file exists and is no longer
    being written to by another program.

    Parameters:
    -----------
    filename: str
        File to check

    Returns:
    --------
    bool
        Indicator for whether file exists or not
    '''
    
    if os.path.isfile(filename):
        size_check_1 = os.path.getsize(filename)
        time.sleep(1)
        size_check_2 = os.path.getsize(filename)
        if size_check_2 > size_check_1:
            return False
        else: 
            return True
    else: 
        return False
    

#------------------------------------------------------------------------------#
#-------------------              Main Program              -------------------#
#------------------------------------------------------------------------------#

def mosaic(dnight, sets, filter, clipFlag):
    '''
    This module creates the mosaic of median filtered images for each data set.

    Parameters:
    -----------
    dnight: str
        Name of data night to process
    sets: list
        List of data sets to process
    filter: str
        Name of photometric filter
    clipFlag: bool
        Whether to use domain clipping (TRUE) or
        rectangle clipping (FALSE) techniques.
    '''
    #set arcpy environment variables part 2/2
    arcpy.CheckOutExtension("Spatial")
    arcpy.env.workspace = filepath.rasters+'scratch_median/'
    arcpy.env.scratchWorkspace = filepath.rasters+'scratch_median'

    #filter paths
    F = {'V':'', 'B':'B/'}
    f = {'V':'', 'B':'b'}

    # Read in the best-fit enxtinction parameters
    extinctionFile = f"{filepath.calibdata}{dnight}/extinction_fit_{filter}.xlsx"
    extinctionData = pd.read_excel(extinctionFile)
    
    for i,s in enumerate(sets):

        # Convert set number to integer
        setnum = int(s[0])

        # Define file paths
        calsetp = f"{filepath.calibdata}{dnight}/S_{setnum:02d}/{F[filter]}"
        gridsetp = f"{filepath.griddata}{dnight}/S_{setnum:02d}/{F[filter]}median/"
        scratchsetp = f"{filepath.rasters}scratch_median/"
        domainsetp = f"{calsetp}/domains/"

        # Remove and/or create gridsetp directory
        if os.path.exists(gridsetp):
            shutil.rmtree(gridsetp, onerror=remove_readonly)
        os.makedirs(gridsetp)

        # Clear scratch directory
        clear_scratch(scratchsetp)
                
        #read in the registered images coordinates
        file = f"{filepath.calibdata}{dnight}/pointerr_{setnum}.txt"
        Obs_AZ, Obs_ALT = n.loadtxt(file, usecols=(3,4)).T
        Obs_AZ[n.where(Obs_AZ>180)] -= 360
        Obs_AZ[35] %= 360
        imnum = len(Obs_AZ)

        # Get zeropoint, platescale, and exposure time for data set
        zeropoint = extinctionData['zeropoint_default'].iloc[setnum-1]
        platescale = extinctionData['avg_scale'].iloc[setnum-1]
        exptime = n.float64(extinctionData['exptime'].iloc[setnum-1])
        
        # Status update
        print(f'{PREFIX}Generating median rasters for {filter}-Band Set {setnum}...')

        #loop through each file in the set
        for w in range(imnum+1):

            v = w+1
            if w == 45:
                w = 35
                Obs_AZ[w] -= 360
            
            # Check that median-filtered TIFF exists
            tiffFile = f'{calsetp}/tiff/median_ib{w+1:03d}.tif'
            while True:
                if check_file(tiffFile): break
                else: time.sleep(1); continue

            arcpy.management.CopyRaster(
                tiffFile, 
                f'ib{v:03d}.tif',
                "DEFAULTS",
                "","","","",
                "16_BIT_UNSIGNED"
            )
            
            #re-define projection to topocentric coordinates
            arcpy.management.DefineProjection(
                f'ib{v:03d}.tif',
                tc(Obs_AZ[w],Obs_ALT[w])
            )
            
            #warp image to remove barrel distortion image
            arcpy.management.Warp(
                f'ib{v:03d}.tif', 
                source_pnt, 
                target_pnt, 
                f'ibw{v:03d}.tif', 
                "POLYORDER3", 
                "BILINEAR"
            )
            set_null_values(f'ibw{v:03d}.tif')

            #reproject into GCS
            arcpy.management.ProjectRaster(
                f'ibw{v:03d}.tif',
                f'fwib{v:03d}.tif',
                geogcs, 
                "BILINEAR", 
                "0.0266"
            )
            set_null_values(f'fwib{v:03d}.tif')
                                       
            # Clip the image
            if clipFlag:
                # Perform domain clipping (clip to true image boundary)

                # Check that clipFile exists first
                clipFile = f'{domainsetp}ib{v:03d}/ib{v:03d}_border'
                while True:
                    if check_file(f"{clipFile}.shp"): break
                    else: time.sleep(1); continue

                arcpy.management.Clip(
                    f"fwib{v:03d}.tif", 
                    "", 
                    f"fcib{v:03d}", 
                    clipFile,
                    "",
                    "ClippingGeometry",
                    "NO_MAINTAIN_EXTENT"
                )
            else:
                # Perform Rectangle clipping
                rectangle = clip_envelope(Obs_AZ, Obs_ALT, w)
                arcpy.management.Clip(
                    f"fwib{v:03d}.tif", rectangle, f"fcib{v:03d}"
                )


            # Status update
            if (v == w+1) & (v % 5 == 0):
                print(f'{PREFIX}{filter}-Band Set {setnum}, {v}/{imnum} rasters complete')
        
        #mosaic raster list must start with an image with max pixel value > 256
        v=1; mstart=1
        while v < (len(Obs_AZ)+1):
            tiff = Image.open(f'{scratchsetp}ib{v:03d}.tif')
            im = n.array(tiff)
            if n.max(im) > 255:
                mstart = v
                break
            v+=1
                        
        #mosaic raster list
        R1 = ';'.join(['fcib%03d' %i for i in range(mstart,47)])
        R2 = ';'.join(['fcib%03d' %i for i in range(1,mstart)])
        R = R1+';'+R2
        
        #mosaic to topocentric coordinate image; save in Griddata\
        print(f"{PREFIX}Mosaicking median rasters for {filter}-Band Set {setnum}...")
        arcpy.management.MosaicToNewRaster(
            R, gridsetp, 'skytopom', geogcs, 
            "32_BIT_FLOAT", "0.0266", "1", "BLEND", "FIRST"
        )

        # Crop lower latitude limit to -6.0 degrees
        arcpy.management.Clip(
            f'{gridsetp}skytopom', 
            "-180.000 -6.000 180.000 90.000", 
            f'{gridsetp}skytopomc'
        )                                      
                                        
        #re-sampling to 0.05 degree resolution
        arcpy.management.Resample(
            gridsetp+'skytopomc',
            gridsetp+'skybright',
            '0.05',
            'BILINEAR'
        )
        
        #convert to magnitudes per square arc second
        print(f"{PREFIX}Converting mosaic to mag/arcsec^2 for {filter}-Band Set {setnum}...")
        psa = 2.5*n.log10((platescale*60)**2) # platescale adjustment
        stm1 = arcpy.sa.Raster(gridsetp + os.sep + 'skybright')
        stm2 = 2.5 * arcpy.sa.Log10(stm1 / exptime)
        skytopomags = zeropoint + psa - stm2
        skytopomags.save(gridsetp+'skybrightmags')
        
        #save mags mosaic to disk
        print(f"{PREFIX}Creating median mosaic layer file for {filter}-Band Set {setnum}...")
        layerName = f"{dnight}_{setnum}_median{f[filter]}"
        layerfile = f"{filepath.griddata}{dnight}/skybrightmags{f[filter]}{setnum}.lyrx"
        symbologyLayer = filepath.rasters+'magnitudes.lyrx'
        arcpy.management.MakeRasterLayer(gridsetp+'skybrightmags', layerName)
        arcpy.management.ApplySymbologyFromLayer(layerName, symbologyLayer)
        arcpy.management.SaveToLayerFile(layerName, layerfile, "ABSOLUTE")

        # Export first dataset's mosaic to JPEG for use in calibreport
        if i == 0:

            # Load in black-background ArcGIS project
            print(f"{PREFIX}Exporting skybright mosaic to JPEG file...")
            blankMap = f"{filepath.maps}blankmap/blankmap.aprx"
            p = arcpy.mp.ArcGISProject(blankMap)

            # Set map scale
            mxd = p.listMaps("Layers")[0]
            lyt = p.listLayouts()[0]
            mf = lyt.listElements()[0]
            mf.camera.scale = 120000000

            # Add sky brightness layer to data frame
            skybrightLayer = arcpy.mp.LayerFile(layerfile)
            mxd.addLayer(skybrightLayer)

            # Save to JPEG from MapView object
            jpegFile = f"{filepath.griddata}{dnight}/S_{setnum:02d}/{F[filter]}data.jpg"
            mv = mxd.defaultView
            mv.exportToJPEG(
                jpegFile,         # output file
                1600,             # width
                1600,             # height
                resolution=96,    # Default is 96
                jpeg_quality=100, # Default is 80, highest quality = 100
            )
        
        #Downscale the raster and save it as a fits file
        file = f"{filepath.griddata}{dnight}/S_{setnum:02d}/{F[filter]}median/skybrightmags"
        arcpy_raster = arcpy.sa.Raster(file)  
        A = arcpy.RasterToNumPyArray(arcpy_raster, "#", "#", "#", -9999)
        A_small = downscale_local_mean(A[:1800,:7200],(25,25)) #72x288
        fname = f"{filepath.griddata}{dnight}/skybrightmags{f[filter]}{setnum}.fits"
        fits.writeto(fname, A_small, overwrite=True)

        # Remove intermediate raster directories & files
        shutil.rmtree(f'{gridsetp}skytopom', onerror=remove_readonly)
        shutil.rmtree(f'{gridsetp}skytopomc', onerror=remove_readonly)
        os.remove(f'{gridsetp}skytopom.aux.xml')
        os.remove(f'{gridsetp}skytopomc.aux.xml')

        # Clear scratch directory again
        clear_scratch(scratchsetp)

        # Final status update
        print(f"{PREFIX}{filter}-Band Set {setnum} median mosaic {pc.CYAN}COMPLETE{pc.END}")
        
    #create mask.tif for horizon masking in the later process
    maskDir = f"{filepath.griddata}{dnight}/mask"
    maskFile = f"{maskDir}/mask.tif"
    if not os.path.exists(maskDir):
        os.makedirs(maskDir)
    if not os.path.isfile(maskFile):
        arcpy.management.CopyRaster(
            f"{gridsetp}skybright",
            maskFile,
            "DEFAULTS",
            "0","0","","",
            "16_BIT_UNSIGNED"
        )

    
if __name__ == "__main__":
    mosaic('FCNA160803', ['1st',],'V')
    pass
