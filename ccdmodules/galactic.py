#-----------------------------------------------------------------------------#
#galactic.py
#
#NPS Night Skies Program
#
#Last updated: 2025/02/05
#
#This script makes the whole sky mosaic of the galactic model according to the 
#time and location of the observed sky. The temporary files generated during the 
#processing stage are stored in filepath.rasters/scratch_galactic folder. The 
#output raster and layer files are stored in file.griddata+dnight.
#
#
#Input: 
#   (1) raster files in the filepath.rasters folder
#   (2) coordinates_%s.txt
#   (3) pointerr_%s.txt
#
#Output:
#   (1) layer files galtopmags%s.lyrx for galactic mosaic
#
#History:
#	Dan Duriscoe -- Created as a module in firstbatchv4vb.py
#	Li-Wei Hung -- Cleaned and improved the code
#   Zach Vanderbosch -- Updated to Python 3.11 and ArcGIS Pro 3.3.1
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
scriptName = 'galactic.py'
PREFIX = f'{pc.GREEN}{scriptName:19s}{pc.END}: '

#-----------------------------------------------------------------------------#
if not os.path.exists(f'{filepath.rasters}scratch_galactic/'):
    os.makedirs(f'{filepath.rasters}scratch_galactic/')

#set arcpy environment variables part 1/2
arcpy.env.rasterStatistics = "STATISTICS 2 2 (-999)"
arcpy.env.overwriteOutput = True
arcpy.env.pyramid = "NONE"

#input rasters for galactic model
galraster  = arcpy.sa.Raster(f'{filepath.rasters}galmagsnew')
galraster1 = arcpy.sa.Raster(f'{filepath.rasters}galmagsnew180')
galrastern = arcpy.sa.Raster(f'{filepath.rasters}galnewnorth')
galrasters = arcpy.sa.Raster(f'{filepath.rasters}galnewsouth')

geogcs = (
    "GEOGCS["
        "'GCS_Sphere_EMEP',"
        "DATUM['D_Sphere_EMEP',"
        "SPHEROID['Sphere_EMEP',6370000.0,0.0]],"
        "PRIMEM['Greenwich',0.0],"
        "UNIT['Degree',0.0174532925199433]"
    "]"
)
    
#-----------------------------------------------------------------------------#
def gal_envelope(lon,lat):
    if abs(lat)<71:
        #expansion angle
        if abs(lat)<55: ang = 22/n.cos(n.deg2rad(lat))
        else: ang = 89
        lon = (lon+90) % 180 - 90
        bond = [lon-ang, lat-19, lon+ang, lat+19] 
    else:
        bond = [-2223549.5, -2223549.5, 2223549.5, 2223549.5]  
    return ' '.join(str(i) for i in bond) #'xmin ymin xmax ymax' 

    
def clip_envelope(AZ, ALT, i):
    if i < 15:
        bond = [AZ[i]-13,-6,AZ[i]+13,ALT[i]+12.9] 
    elif i < 30: 
        bond = [AZ[i]-13,ALT[i]-12.7,AZ[i]+13,ALT[i]+12.6] 
    elif i < 40:
        bond = [AZ[i]-18.6,ALT[i]-12.7,AZ[i]+18.6,ALT[i]+12.6] 
    elif i < 45:
        bond = [AZ[i]-39.6,ALT[i]-12.7,AZ[i]+39.6,ALT[i]+12.7] 
    return ' '.join(str(i) for i in bond) #'xmin ymin xmax ymax' 


def tc(lon,lat):
    '''
    Returns the topocentric coordinate setting in WKT format
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
    
    
def get_galgn(lon,lat):
    rectangle = gal_envelope(lon,lat)
    if abs(lat)<71:
        if abs(lon)<90: arcpy.management.Clip(galraster,rectangle,"galclip.tif")
        else: arcpy.management.Clip(galraster1, rectangle, "galclip.tif")
        lon = (lon+90) % 180 - 90
        p = [tc(lon,lat),"BILINEAR", "6000"]
        arcpy.management.ProjectRaster("galclip.tif", "galgn.tif", *p)
    else: 
        p = [tc(lon,lat),"BILINEAR", "6000"]
        if lat>=71:
            arcpy.management.ProjectRaster(galrastern, "galtemp.tif", *p)
        else:
            arcpy.management.ProjectRaster(galrasters, "galtemp.tif", *p)
        arcpy.management.Clip('galtemp.tif', rectangle, 'galgn.tif')

  
def remove_readonly(func, path, excinfo):
    '''
    Error-catching function to handle removal of read-only folders
    '''
    os.chmod(path, stat.S_IWRITE)
    func(path)
  
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


def check_file(filename):
    ''' 
    Check that a file exists and is no longer
    being written to by another program.
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
 

def mosaic(dnight, sets):
    '''
    This module creates the mosaic of the galactic model for each data set.
    '''
    # Set arcpy environment variables part 2/2
    arcpy.env.workspace = f'{filepath.rasters}scratch_galactic/'
    arcpy.env.scratchWorkspace = f'{filepath.rasters}scratch_galactic'

    for s in sets:

        # Define file paths
        calsetp = f"{filepath.calibdata}{dnight}/"
        gridsetp = f"{filepath.griddata}{dnight}/S_0{s[0]}/gal/"
        scratchsetp = f"{filepath.rasters}scratch_galactic/"
        domainsetp = f"{calsetp}/S_0{s[0]}/domains/"

        # Remove and/or create gridsetp directory
        if os.path.exists(gridsetp):
            shutil.rmtree(gridsetp, onerror=remove_readonly)
        os.makedirs(gridsetp)

        # Clear scratch directory
        clear_scratch(scratchsetp)
        
        # Read in the galactic coordinates from coordinates_%s.txt
        file = f'{calsetp}coordinates_{s[0]}.txt'
        Gal_ang, Gal_l, Gal_b = n.loadtxt(file,usecols=(1,2,3)).T
        
        # Read in the registered images coordinates
        file = f'{calsetp}pointerr_{s[0]}.txt'
        Obs_AZ, Obs_ALT = n.loadtxt(file, usecols=(3,4)).T
        Obs_AZ[n.where(Obs_AZ>180)] -= 360
        Obs_AZ[35] %= 360
        imnum = len(Obs_AZ)
        
        # Status update
        print(f'{PREFIX}Generating galactic rasters for Set {s[0]}...')

        # Loop through each file in the set
        for w in range(imnum+1):
            
            v = w+1
            if w == 45:
                w = 35
                Obs_AZ[w] -= 360
                
            get_galgn(Gal_l[w], Gal_b[w])
        
            # Rotate by galactic angle
            arcpy.management.Rotate(
                'galgn.tif', 
                'rotaterasterg.tif', 
                str(Gal_ang[w]), 
                "0 0",
                "BILINEAR"
            )
                                    
            # Re-define projection to topocentric coordinates
            arcpy.management.DefineProjection(
                'rotaterasterg.tif',
                tc(Obs_AZ[w],Obs_ALT[w])
            )
                                            
            # Reproject into GCS
            arcpy.management.ProjectRaster(
                'rotaterasterg.tif', 
                f'gal{v:02d}.tif', 
                geogcs,
                "BILINEAR",
                "0.05"
            )
                                        
            # Clip to image boundary
            # rectangle = clip_envelope(Obs_AZ, Obs_ALT, w)
            # arcpy.management.Clip(f"gal{v:02d}.tif", rectangle, f"gali{v:02d}")

            # Check that clipFile exists first
            clipFile = f'{domainsetp}ib{v:03d}/ib{v:03d}_border'
            while True:
                if check_file(f"{clipFile}.shp"): break
                else: time.sleep(1); continue

            arcpy.management.Clip(
                f"gal{v:02d}.tif", 
                "", 
                f"gali{v:02d}", 
                clipFile,
                "0",
                "ClippingGeometry",
                "NO_MAINTAIN_EXTENT"
            )

            # Status update
            if (v == w+1) & (v % 5 == 0):
                print(f'{PREFIX}Set {s[0]}, {v}/{imnum} rasters complete')
        
        # Mosaic to topocentric coordinate model; save in Griddata\
        print(f"{PREFIX}Mosaicking galactic rasters for Set {s[0]}...")
        R = ';'.join([f'gali{i:02d}' for i in range(1,47)])
        arcpy.management.MosaicToNewRaster(
            R, gridsetp, 'galtopmagsuc', geogcs, 
            "32_BIT_FLOAT", "0.05", "1", "BLEND", "FIRST"
        )

        # Crop lower latitude limit to -6.0 degrees
        arcpy.management.Clip(
            f'{gridsetp}galtopmagsuc', 
            "-180.000 -6.000 180.000 90.000", 
            f'{gridsetp}galtopmags'
        )

        # Create Raster layer, add magnitudes symbology, and save layer to file
        print(f"{PREFIX}Creating galactic mosaic layer file for Set {s[0]}...")
        layerfile = f'{filepath.griddata}{dnight}/galtopmags{s[0]}.lyrx'
        symbologyFile = f'{filepath.rasters}magnitudes.lyrx'
        arcpy.management.MakeRasterLayer(gridsetp+'galtopmags', 'galtoplyr')
        arcpy.management.ApplySymbologyFromLayer('galtoplyr', symbologyFile)
        arcpy.management.SaveToLayerFile('galtoplyr', layerfile, "ABSOLUTE")
        
        # Downscale the raster and save it as a fits file
        file = f'{gridsetp}galtopmags'
        arcpy_raster = arcpy.sa.Raster(file)
        A = arcpy.RasterToNumPyArray(arcpy_raster, "#", "#", "#", -9999)
        A_small = downscale_local_mean(A[:1800,:7200],(25,25)) #72x288
        fname = f'{filepath.griddata}{dnight}/galtopmags{s[0]}.fits'
        fits.writeto(fname, A_small, overwrite=True)

        # Remove intermediate raster directories & files
        shutil.rmtree(f'{gridsetp}galtopmagsuc', onerror=remove_readonly)
        os.remove(f'{gridsetp}galtopmagsuc.aux.xml')

        # Clear scratch directory again
        clear_scratch(scratchsetp)

        # Status update
        print(f"{PREFIX}Set {s[0]} galactic mosaic {pc.CYAN}COMPLETE{pc.END}")

if __name__ == "__main__":
    mosaic('FCNA160803', ['1st',])
    pass








