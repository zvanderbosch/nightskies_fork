#-----------------------------------------------------------------------------#
#zodiacal.py
#
#NPS Night Skies Program
#
#Last updated: 2025/07/24
#
#This script makes the whole sky mosaic of the zodiacal model according to the 
#time and location of the observed sky. The temporary files generated during the 
#processing stage are stored in filepath.rasters/scratch_zodiacal folder. The 
#output raster and layer files are stored in file.griddata+dnight.
#
#
#Input: 
#   (1) pointerr_<DATASET>.txt
#           Pointing error data for image Alt/Az coordinates
#           (filepath.calibdata/DATANIGHT)
#   (2) coordinates_<DATASET>.txt
#           Solved image coordinates for Ecliptic coordinates
#           (filepath.calibdata/DATANIGHT)
#   (3) zodiacal_01
#           Raster dataset for zodiacal light model (eastern sky)
#           (filepath.rasters)
#   (4) zodiacal_180
#           Raster dataset for zodiacal light model (western sky)
#           (filepath.rasters)
#
#Output:
#   (1) zodtopmags
#           Zodiacal light model mosaic
#           (filepath.griddata/DATANIGHT/S_0#/zod)
#   (2) zodtopmags<DATASET>.lyrx 
#           Layer file for zodiacal light model mosaic
#           (filepath.griddata/DATANIGHT)
#
#History:
#	Dan Duriscoe -- Created as a module in firstbatchv4vb.py
#	Li-Wei Hung -- Cleaned and improved the code
#   Zach Vanderbosch -- Updated to Python 3.11 and ArcGIS Pro 3.5.2
#
#-----------------------------------------------------------------------------#
from glob import glob
from astropy.io import fits
from skimage.transform import downscale_local_mean

import arcpy
import numpy as n
import os
import time
import stat
import shutil

# Local Source
import filepath
import printcolors as pc

# Print status prefix
scriptName = 'zodiacal.py'
PREFIX = f'{pc.GREEN}{scriptName:19s}{pc.END}: '

#-----------------------------------------------------------------------------#
if not os.path.exists(filepath.rasters+'scratch_zodiacal/'):
    os.makedirs(filepath.rasters+'scratch_zodiacal/')

#set arcpy environment variables part 1/2
arcpy.env.rasterStatistics = "NONE"
arcpy.env.overwriteOutput = True
arcpy.env.pyramid = "NONE"

#input rasters for zodiacal model
zodraster = arcpy.sa.Raster(filepath.rasters+'zodiacal_01')
zodraster1 = arcpy.sa.Raster(filepath.rasters+'zodiacal_180')
    
geogcs = (
    "GEOGCS["
        "'GCS_Sphere_EMEP',"
        "DATUM['D_Sphere_EMEP',"
        "SPHEROID['Sphere_EMEP',6370000.0,0.0]],"
        "PRIMEM['Greenwich',0.0],"
        "UNIT['Degree',0.0174532925199433]"
    "]"
)
          
#------------------------------------------------------------------------------#
#-------------------            Define Functions            -------------------#
#------------------------------------------------------------------------------#

def zod_envelope(lon,lat):
    '''
    Returns clipping boundaries for the zodiacal model

    Parameters:
    -----------
    lon: float
        Longitude coordinate
    lat: float
        Latitude coordinate

    Returns:
    --------
    boundary: str
        Coordinate boundaries for galactic model
    '''
    if abs(lat)<71:
        #expansion angle
        if abs(lat)<55: ang = 22/n.cos(n.deg2rad(lat))
        else: ang = 89
        lon = (lon+90) % 180 - 90
        bond = [lon-ang, lat-19, lon+ang, lat+19] 
    else:
        bond = [-2000000, -2000000, 2000000, 2000000]  
    boundary = ' '.join(str(i) for i in bond) #'xmin ymin xmax ymax'
    return boundary 
    
    
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
    

def get_zodgn(lon,lat):
    '''
    Generates zodiacal model projected and clipped 
    to local site.

    Parameters:
    -----------
    lon: float
        Longitude coordinate
    lat: float
        Latitude coordinate
    '''
    rectangle = zod_envelope(lon,lat)
    if abs(lat)<71:
        if abs(lon)<90: arcpy.management.Clip(zodraster,rectangle,"zodclip.tif")
        else: arcpy.management.Clip(zodraster1, rectangle, "zodclip.tif")
        lon = (lon+90) % 180 - 90
        p = [tc(lon,lat),"BILINEAR", "6000"]
        arcpy.management.ProjectRaster("zodclip.tif", "zodgn.tif", *p)
    else: 
        p = [tc(lon,lat),"BILINEAR", "6000"]
        arcpy.management.ProjectRaster(zodraster, "zodtemp.tif", *p)
        arcpy.management.Clip('zodtemp.tif', rectangle, 'zodgn.tif')


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

def mosaic(dnight, sets):
    '''
    This module creates the mosaic of the zodiacal model for each data set.
    
    Parameters:
    -----------
    dnight: str
        Name of data night to process
    sets: list
        List of data sets to process
    '''
    #set arcpy environment variables part 2/2
    arcpy.CheckOutExtension("Spatial")
    arcpy.env.workspace = filepath.rasters+'scratch_zodiacal/'
    arcpy.env.scratchWorkspace = filepath.rasters+'scratch_zodiacal'

    for s in sets:

        # Conver set number to integer
        setnum = int(s[0])

        # Define file paths
        calsetp = f"{filepath.calibdata}{dnight}/S_{setnum:02d}/"
        gridsetp = f"{filepath.griddata}{dnight}/S_{setnum:02d}/zod/"
        scratchsetp = f"{filepath.rasters}scratch_zodiacal/"
        domainsetp = f"{calsetp}domains/"

        if os.path.exists(gridsetp):
            shutil.rmtree(gridsetp, onerror=remove_readonly)
        os.makedirs(gridsetp)

        # Clear out scratch directory
        clear_scratch(scratchsetp)
        
        #read in the zodiacal coordinates from coordinates_%s.txt
        file = f"{filepath.calibdata}{dnight}/coordinates_{setnum}.txt"
        Ecl_ang, Ecl_l, Ecl_b = n.loadtxt(file,usecols=(4,5,6)).T
        
        #read in the registered images coordinates
        file = f"{filepath.calibdata}{dnight}/pointerr_{setnum}.txt"
        Obs_AZ, Obs_ALT = n.loadtxt(file, usecols=(3,4)).T
        Obs_AZ[n.where(Obs_AZ>180)] -= 360
        Obs_AZ[35] %= 360
        imnum = len(Obs_AZ)
        
        # Status update
        print(f'{PREFIX}Generating zodiacal rasters for Set {setnum}...')

        #loop through each file in the set
        for w in range(imnum+1):

            v = w+1
            if w == 45:
                w = 35
                Obs_AZ[w] -= 360
            
            get_zodgn(Ecl_l[w], Ecl_b[w])
            
            #rotate by zodiacal angle
            arcpy.management.Rotate(
                'zodgn.tif', 
                'rotateraster.tif', 
                str(Ecl_ang[w]), 
                "0 0",
                "BILINEAR"
            )
                                    
            #re-define projection to topocentric coordinates
            arcpy.management.DefineProjection(
                'rotateraster.tif',
                tc(Obs_AZ[w],Obs_ALT[w])
            )

            #reproject into GCS
            arcpy.management.ProjectRaster(
                'rotateraster.tif', 
                'zod%02d.tif'%v, 
                geogcs,
                "BILINEAR",
                "0.1"
            )

            #clip to image boundary
            # rectangle = clip_envelope(Obs_AZ, Obs_ALT, w)
            # arcpy.management.Clip("zod%02d.tif"%v, rectangle, "zodi%02d"%v)

            # Check that clipFile exists
            clipFile = f'{domainsetp}ib{v:03d}/ib{v:03d}_border'
            while True:
                if check_file(f"{clipFile}.shp"): break
                else: time.sleep(1); continue

            arcpy.management.Clip(
                f"zod{v:02d}.tif", 
                "", 
                f"zodi{v:02d}", 
                clipFile,
                "0",
                "ClippingGeometry",
                "NO_MAINTAIN_EXTENT"
            )

            # Status update
            if (v == w+1) & (v % 5 == 0):
                print(f'{PREFIX}Set {setnum}, {v}/{imnum} rasters complete')
            
        #Mosaic to topocentric coordinate model; save in Griddata\
        print(f"{PREFIX}Mosaicking zodiacal rasters for Set {setnum}...")
        R = ';'.join(['zodi%02d' %i for i in range(1,47)])
        arcpy.management.MosaicToNewRaster(
            R, gridsetp, 'zodtopo', geogcs, 
            "32_BIT_FLOAT", "0.1", "1", "BLEND", "FIRST"
        )

        # Crop lower latitude limit to -6.0 degrees
        arcpy.management.Clip(
            f'{gridsetp}zodtopo', 
            "-180.000 -6.000 180.000 90.000", 
            f'{gridsetp}zodtopoc'
        )
                                        
        #re-sampling to 0.05 degree resolution
        gridname = gridsetp + "zodtopmags"
        arcpy.management.Resample(
            gridsetp+'zodtopoc',
            gridname,'0.05',
            'BILINEAR'
        )
    
        #Create Raster layer, add magnitudes symbology, and save layer to file
        print(f"{PREFIX}Creating zodiacal mosaic layer file for Set {setnum}...")
        layerfile = f"{filepath.griddata}{dnight}/zodtopmags{setnum}.lyrx"
        symbologyLayer = f"{filepath.rasters}magnitudes.lyrx"
        arcpy.management.MakeRasterLayer(gridsetp+'zodtopmags', 'zodtoplyr')
        arcpy.management.ApplySymbologyFromLayer('zodtoplyr', symbologyLayer)
        arcpy.management.SaveToLayerFile('zodtoplyr', layerfile, "ABSOLUTE")
        
        #Downscale the raster and save it as a fits file
        file = f"{filepath.griddata}{dnight}/S_{setnum:02d}/zod/zodtopmags"
        arcpy_raster = arcpy.sa.Raster(file)  
        A = arcpy.RasterToNumPyArray(arcpy_raster, "#", "#", "#", -9999)
        A_small = downscale_local_mean(A[:1800,:7200],(25,25)) #72x288
        fname = f"{filepath.griddata}{dnight}/zodtopmags{setnum}.fits"
        fits.writeto(fname, A_small, overwrite=True)

        # Remove intermediate raster directories & files
        shutil.rmtree(f'{gridsetp}zodtopo', onerror=remove_readonly)
        shutil.rmtree(f'{gridsetp}zodtopoc', onerror=remove_readonly)
        os.remove(f'{gridsetp}zodtopo.aux.xml')
        os.remove(f'{gridsetp}zodtopoc.aux.xml')

        # Clear scratch directory again
        clear_scratch(scratchsetp)

        # Final status update
        print(
            f"{PREFIX}Set {setnum} zodiacal mosaic {pc.CYAN}COMPLETE{pc.END}"
        )

    
if __name__ == "__main__":
    mosaic('FCNA160803', ['1st',])
    pass
