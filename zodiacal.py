#-----------------------------------------------------------------------------#
#zodiacal.py
#
#NPS Night Skies Program
#
#Last updated: 2025/02/05
#
#This script makes the whole sky mosaic of the zodiacal model according to the 
#time and location of the observed sky. The temporary files generated during the 
#processing stage are stored in filepath.rasters/scratch_zodiacal folder. The 
#output raster and layer files are stored in file.griddata+dnight.
#
#
#Input: 
#   (1) raster files in the filepath.rasters folder
#   (2) coordinates_%s.txt
#   (3) pointerr_%s.txt
#
#Output:
#   (1) layer files zodtopmags%s.lyrx for zodiacal mosaic
#
#History:
#	Dan Duriscoe -- Created as a module in firstbatchv4vb.py
#	Li-Wei Hung -- Cleaned and improved the code
#   Zach Vanderbosch -- Updated to Python 3.11 and ArcGIS Pro 3.3.1
#
#-----------------------------------------------------------------------------#
from astropy.io import fits
from tqdm import trange
from skimage.transform import downscale_local_mean

import arcpy
import numpy as n
import os
import stat
import shutil

# Local Source
import filepath  

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
          
#-----------------------------------------------------------------------------#
def zod_envelope(lon,lat):
    if abs(lat)<71:
        #expansion angle
        if abs(lat)<55: ang = 22/n.cos(n.deg2rad(lat))
        else: ang = 89
        lon = (lon+90) % 180 - 90
        bond = [lon-ang, lat-19, lon+ang, lat+19] 
    else:
        bond = [-2000000, -2000000, 2000000, 2000000]  
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
    

def get_zodgn(lon,lat):
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
    This module creates the mosaic of the zodiacal model for each data set.
    '''
    #set arcpy environment variables part 2/2
    arcpy.CheckOutExtension("Spatial")
    arcpy.env.workspace = filepath.rasters+'scratch_zodiacal/'
    arcpy.env.scratchWorkspace = filepath.rasters+'scratch_zodiacal'

    for s in sets:

        # Define file paths
        calsetp = filepath.calibdata+dnight+'/S_0%s/' %s[0]
        gridsetp = filepath.griddata+dnight+'/S_0%s/zod/' %s[0]
        scratchsetp = f"{filepath.rasters}scratch_zodiacal/"
        domainsetp = f"{calsetp}domains/"

        if os.path.exists(gridsetp):
            shutil.rmtree(gridsetp, onerror=remove_readonly)
        os.makedirs(gridsetp)

        # Clear out scratch directory
        clear_scratch(scratchsetp)
        
        #read in the zodiacal coordinates from coordinates_%s.txt
        file = filepath.calibdata+dnight+'/coordinates_%s.txt'%s[0]
        Ecl_ang, Ecl_l, Ecl_b = n.loadtxt(file,usecols=(4,5,6)).T
        
        #read in the registered images coordinates
        file = filepath.calibdata+dnight+'/pointerr_%s.txt' %s[0]
        Obs_AZ, Obs_ALT = n.loadtxt(file, usecols=(3,4)).T
        Obs_AZ[n.where(Obs_AZ>180)] -= 360
        Obs_AZ[35] %= 360
        
        #loop through each file in the set
        print(f'Generating zodiacal images for Set {s[0]}...')
        for w in trange(len(Obs_AZ)+1):
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
            clipFile = f'{domainsetp}ib{v:03d}/ib{v:03d}_border'
            arcpy.management.Clip(
                f"zod{v:02d}.tif", 
                "", 
                f"zodi{v:02d}", 
                clipFile,
                "0",
                "ClippingGeometry",
                "NO_MAINTAIN_EXTENT"
            )
            
        #Mosaic to topocentric coordinate model; save in Griddata\
        print("Mosaicking into all sky zodiacal model...")
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
        print("Creating layer files for zodiacal mosaic...")
        layerfile = filepath.griddata+dnight+'/zodtopmags%s.lyrx' %s[0]
        symbologyLayer = filepath.rasters+'magnitudes.lyrx'
        arcpy.management.MakeRasterLayer(gridsetp+'zodtopmags', 'zodtoplyr')
        arcpy.management.ApplySymbologyFromLayer('zodtoplyr', symbologyLayer)
        arcpy.management.SaveToLayerFile('zodtoplyr', layerfile, "ABSOLUTE")
        
        #Downscale the raster and save it as a fits file
        file = filepath.griddata+dnight+"/S_0"+s[0]+"/zod/zodtopmags"
        arcpy_raster = arcpy.sa.Raster(file)  
        A = arcpy.RasterToNumPyArray(arcpy_raster, "#", "#", "#", -9999)
        A_small = downscale_local_mean(A[:1800,:7200],(25,25)) #72x288
        fname = filepath.griddata+dnight+'/zodtopmags%s.fits' %s[0]
        fits.writeto(fname, A_small, overwrite=True)

    
if __name__ == "__main__":
    mosaic('FCNA160803', ['1st',])
    pass
