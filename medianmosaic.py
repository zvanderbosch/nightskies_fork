#-----------------------------------------------------------------------------#
#medianmosaic.py
#
#NPS Night Skies Program
#
#Last updated: 2017/02/02
#
#This script makes the whole sky mosaic from the median filtered images 
#according to the location of observed sky. The temporary files generated during
#the processing stage are stored in filepath.rasters/scratch_median folder. The 
#output raster and layer files are stored in file.griddata+dnight.
#
#
#Input: 
#   (1) median filtered tiff files in the filepath.calibdata
#   (2) pointerr_%s.txt
#   (3) extinction_fit_%s.txt
#   (4) raster files in the filepath.rasters folder
#
#Output:
#   (1) layer files skybrightmags%s.lyr for median mosaic
#   (2) mask.tif for making the horizontal mask in the later process
#
#History:
#	Dan Duriscoe -- Created as a module in firstbatchv4vb.py
#	Li-Wei Hung -- Cleaned and improved the code
#
#-----------------------------------------------------------------------------#
from astropy.io import fits
from tqdm import trange
from PIL import Image
from skimage.transform import downscale_local_mean

import arcpy
import numpy as n
import os
import stat
import shutil

# Local Source
import filepath  

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
          
#-----------------------------------------------------------------------------#
def clip_envelope(AZ, ALT, i):
    '''
    Calculate the four coordinates defining the minimum bounding rectangle to be clipped
    by the arcpy.Clip_management() function. 

    Returns coordinates are defined in this order: 
    X-Minimum, Y-Minimum, X-Maximum, Y-Maximum.
    '''
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


def set_null_values(rasterFile):
    '''
    Function to set values within raster File to NoData
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
    

def mosaic(dnight, sets, filter):
    '''
    This module creates the mosaic of median filtered images for each data set.
    '''
    #set arcpy environment variables part 2/2
    arcpy.CheckOutExtension("Spatial")
    arcpy.env.workspace = filepath.rasters+'scratch_median/'
    arcpy.env.scratchWorkspace = filepath.rasters+'scratch_median'

    #filter paths
    F = {'V':'', 'B':'B/'}
    f = {'V':'', 'B':'b'}
    
    for s in sets:

        # Define file paths
        calsetp = f"{filepath.calibdata}{dnight}/S_0{s[0]}/{F[filter]}"
        gridsetp = f"{filepath.griddata}{dnight}/S_0{s[0]}/{F[filter]}median/"
        scratchsetp = f"{filepath.rasters}scratch_median/"
        domainsetp = f"{calsetp}/domains/"

        # Remove and/or create gridsetp directory
        if os.path.exists(gridsetp):
            shutil.rmtree(gridsetp, onerror=remove_readonly)
        os.makedirs(gridsetp)

        # Clear scratch directory
        clear_scratch(scratchsetp)
                
        #read in the registered images coordinates
        file = filepath.calibdata+dnight+'/pointerr_%s.txt' %s[0]
        Obs_AZ, Obs_ALT = n.loadtxt(file, usecols=(3,4)).T
        Obs_AZ[n.where(Obs_AZ>180)] -= 360
        Obs_AZ[35] %= 360
        
        #read in the best-fit zeropoint and plate scale
        file = filepath.calibdata+dnight+'/extinction_fit_%s.txt' %filter
        zeropoint, platescale, exptime = n.loadtxt(file, usecols=(2,8,9), unpack=True, ndmin=2)
        
        #loop through each file in the set
        print(f'Generating median images for Set {s[0]}...')
        for w in trange(len(Obs_AZ)+1):

            v = w+1
            if w == 45:
                w = 35
                Obs_AZ[w] -= 360
            
            arcpy.management.CopyRaster(
                f'{calsetp}/tiff/median_ib{w+1:03d}.tif', 
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
                                       
            # Clip to image boundary
            # rectangle = clip_envelope(Obs_AZ, Obs_ALT, w)
            # arcpy.management.Clip(
            #     f"fwib{v:03d}.tif", 
            #     rectangle, 
            #     f"fcib{v:03d}"
            # )
            clipFile = f'{domainsetp}ib{v:03d}/ib{v:03d}_border'
            arcpy.management.Clip(
                f"fwib{v:03d}.tif", 
                "", 
                f"fcib{v:03d}", 
                clipFile,
                "",
                "ClippingGeometry",
                "NO_MAINTAIN_EXTENT"
            )
        
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
        print("Mosaicking into all sky median image...")
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
        print("Converting the mosaic to mag per squard arcsec...")
        psa = 2.5*n.log10((platescale[int(s[0])-1]*60)**2) # platescale adjustment
        stm1 = arcpy.sa.Raster(gridsetp + os.sep + 'skybright')
        stm2 = stm1 / exptime[0]
        stm3 = 2.5 * arcpy.sa.Log10(stm2)
        skytopomags = zeropoint[int(s[0])-1] + psa - stm3
        
        #save mags mosaic to disk
        skytopomags.save(gridsetp+'skybrightmags')
    
        print("Creating layer files for median mosaic...")
        layerName = dnight+'_%s_median%s'%(s[0],f[filter])
        layerfile = filepath.griddata+dnight+'/skybrightmags%s%s.lyr'%(f[filter],s[0])
        symbologyLayer = filepath.rasters+'magnitudes.lyr'
        arcpy.management.MakeRasterLayer(gridsetp+'skybrightmags', layerName)
        arcpy.management.ApplySymbologyFromLayer(layerName, symbologyLayer)
        arcpy.management.SaveToLayerFile(layerName, layerfile, "ABSOLUTE")
        
        #Downscale the raster and save it as a fits file
        file = filepath.griddata+dnight+'/S_0%s/%smedian/skybrightmags' %(s[0],F[filter])
        arcpy_raster = arcpy.sa.Raster(file)  
        A = arcpy.RasterToNumPyArray(arcpy_raster, "#", "#", "#", -9999)
        A_small = downscale_local_mean(A[:1800,:7200],(25,25)) #72x288
        fname = filepath.griddata+dnight+'/skybrightmags%s%s.fits'%(f[filter],s[0])
        fits.writeto(fname, A_small, overwrite=True)
        
    #create mask.tif for horizon masking in the later process
    mask = filepath.griddata+dnight+'/mask.tif'
    if not os.path.isfile(mask):
        arcpy.management.CopyRaster(
            gridsetp+'skybright',
            mask,
            "DEFAULTS",
            "0","0","","",
            "16_BIT_UNSIGNED"
        )

    
if __name__ == "__main__":
    mosaic('FCNA160803', ['1st',],'V')
    pass
