#-----------------------------------------------------------------------------#
#galactic.py
#
#NPS Night Skies Program
#
#Last updated: 2016/12/23
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
#   (1) layer files galtopmags%s.lyr for galactic mosaic
#
#History:
#	Dan Duriscoe -- Created as a module in firstbatchv4vb.py
#	Li-Wei Hung -- Cleaned and improved the code
#   Zach Vanderbosch -- Updated to Python3 and ArcGIS Pro
#
#-----------------------------------------------------------------------------#
from astropy.io import fits
from skimage.transform import downscale_local_mean
from tqdm import trange

import arcpy
import numpy as n
import os
import stat
import shutil

# Local Source
import filepath  

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

geogcs = "GEOGCS['GCS_Sphere_EMEP',\
          DATUM['D_Sphere_EMEP',SPHEROID['Sphere_EMEP',6370000.0,0.0]],\
          PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]]"
    
#-----------------------------------------------------------------------------#
def gal_envelope(lon,lat):
    if abs(lat)<71:
        #expension angle
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
    '''Returns the topocentric coordinate setting'''
    return "PROJCS['gnomonic',%s,PROJECTION['Gnomonic'],\
    PARAMETER['False_Easting',0.0],PARAMETER['False_Northing',0.0],\
    PARAMETER['Longitude_Of_Center',%s],PARAMETER['Latitude_Of_Center',%s],\
    UNIT['Meter',1.0]]"%(geogcs,str(lon),str(lat))
    
    
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
        

def mosaic(dnight, sets):
    '''
    This module creates the mosaic of the galactic model for each data set.
    '''
    # Set arcpy environment variables part 2/2
    arcpy.env.workspace = f'{filepath.rasters}scratch_galactic/'
    arcpy.env.scratchWorkspace = f'{filepath.rasters}scratch_galactic'

    for s in sets:

        # Clear out scratch directory
        clear_scratch(f'{filepath.rasters}scratch_galactic/')

        # File paths
        calsetp = f"{filepath.calibdata}{dnight}/"
        gridsetp = f"{filepath.griddata}{dnight}/S_0{s[0]}/gal/"
        if os.path.exists(gridsetp):
            shutil.rmtree(gridsetp, onerror=remove_readonly)
        os.makedirs(gridsetp)
        
        # Read in the galactic coordinates from coordinates_%s.txt
        file = f'{calsetp}coordinates_{s[0]}.txt'
        Gal_ang, Gal_l, Gal_b = n.loadtxt(file,usecols=(1,2,3)).T
        
        # Read in the registered images coordinates
        file = f'{calsetp}pointerr_{s[0]}.txt'
        Obs_AZ, Obs_ALT = n.loadtxt(file, usecols=(3,4)).T
        Obs_AZ[n.where(Obs_AZ>180)] -= 360
        Obs_AZ[35] %= 360
        
        # Loop through each file in the set
        print('Generating galactic images for Set {s[0]}...')
        for w in trange(len(Obs_AZ)+1):
            v = w+1
            if w == 45:
                w = 35
                Obs_AZ[w] -= 360
                
            get_galgn(Gal_l[w], Gal_b[w])
        
            # Rotate by galctic angle
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
            rectangle = clip_envelope(Obs_AZ, Obs_ALT, w)
            arcpy.management.Clip(f"gal{v:02d}.tif", rectangle, f"gali{v:02d}")
        
        # Mosaic to topocentric coordinate model; save in Griddata\
        print("Mosaicking into all sky galactic model...")
        R = ';'.join([f'gali{i:02d}' for i in range(1,47)])
        arcpy.management.MosaicToNewRaster(
            R, gridsetp, 'galtopmags', geogcs, 
            "32_BIT_FLOAT", "0.05", "1", "BLEND", "FIRST"
        )

        # Create Raster layer, add magnitudes symbology, and save layer to file
        print("Creating layer files for galactic mosaic...")
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

if __name__ == "__main__":
    mosaic('FCNA160803', ['1st',])
    pass








