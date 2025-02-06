#-----------------------------------------------------------------------------#
#fullmosaic.py
#
#NPS Night Skies Program
#
#Last updated: 2025/02/05
#
#This script makes the whole sky mosaic from the full resolution images 
#according to the location of observed sky. The temporary files generated during
#the processing stage are stored in filepath.rasters/scratch_fullres folder. The 
#output raster and layer files are stored in file.griddata+dnight.
#
#
#Input: 
#   (1) full resolution tiff files in the filepath.calibdata
#   (2) pointerr_%s.txt
#   (3) extinction_fit_%s.txt
#   (4) raster files in the filepath.rasters folder
#
#Output:
#   (1) layer files skytopomags%s.lyrx for full-resolution mosaic
#
#History:
#	Dan Duriscoe -- Created as a module in firstbatchv4vb.py
#	Li-Wei Hung -- Cleaned and improved the code
#   Zach Vanderbosch -- Updated to Python 3.11 and ArcGIS Pro 3.3.1
#
#-----------------------------------------------------------------------------#

from PIL import Image

import arcpy
import numpy as n
import os
import stat
import shutil

# Local Source
import filepath
import printcolors as pc

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
          
#-----------------------------------------------------------------------------#
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
  
  
def clear_dir(dir_path):
    '''
    Function to clear out all files and folders from
    the specified directory.
    '''
    for root, dirs, files in os.walk(dir_path, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.chmod(os.path.join(root, name), stat.S_IWRITE)
            os.rmdir(os.path.join(root, name))
    

def mosaic(dnight, sets, filter):
    '''
    This module creates the mosaic of full-resolution images for each data set.
    '''
    #set arcpy environment variables part 2/2
    arcpy.CheckOutExtension("Spatial")
    arcpy.env.workspace = filepath.rasters+'scratch_fullres/'
    arcpy.env.scratchWorkspace = filepath.rasters+'scratch_fullres'

    #filter paths
    F = {'V':'', 'B':'B/'}
    f = {'V':'', 'B':'b'}
    
    for s in sets:

        # Define file paths
        calsetp = f"{filepath.calibdata}{dnight}/S_0{s[0]}/{F[filter]}"
        gridsetp = f"{filepath.griddata}{dnight}/S_0{s[0]}/{F[filter]}fullres/"
        scratchsetp = f"{filepath.rasters}scratch_fullres/"
        domainsetp = f"{calsetp}/domains/"

        # Remove and/or create gridsetp directory
        if os.path.exists(gridsetp):
            shutil.rmtree(gridsetp, onerror=remove_readonly)
        os.makedirs(gridsetp)

        # Create domainsetp if non-existent, else clear out directory
        if not os.path.exists(domainsetp):
            os.makedirs(domainsetp)
        else:
            clear_dir(domainsetp)

        # Clear out the scratch directory
        clear_dir(scratchsetp)
                
        #read in the registered images coordinates
        file = filepath.calibdata+dnight+'/pointerr_%s.txt' %s[0]
        Obs_AZ, Obs_ALT = n.loadtxt(file, usecols=(3,4)).T
        Obs_AZ[n.where(Obs_AZ>180)] -= 360
        Obs_AZ[35] %= 360
        imnum = len(Obs_AZ)
        
        #read in the best-fit zeropoint and plate scale
        file = filepath.calibdata+dnight+'/extinction_fit_%s.txt' %filter
        zeropoint, platescale, exptime = n.loadtxt(
            file, usecols=(2,8,9), unpack=True, ndmin=2
        )
        
        # Status Update
        print(
            f'{pc.GREEN}fullmosaic.py  {pc.END}'
            f': Generating fullres rasters for {filter}-Band Set {s[0]}...'
        )

        #loop through each file in the set
        for w in range(imnum+1):

            v = w+1
            if w == 45:
                w = 35
                Obs_AZ[w] -= 360

            # Copy TIFF file to scratch directory
            arcpy.management.CopyRaster(
                f'{calsetp}/tiff/ib{w+1:03d}.tif', 
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
                print(
                    f'{pc.GREEN}fullmosaic.py  {pc.END}'
                    f': {filter}-Band Set {s[0]}, {v}/{imnum} rasters complete'
                )
            
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
        print(
            f"{pc.GREEN}fullmosaic.py  {pc.END}"
            f": Mosaicking into all sky full-resolution image for {filter}-Band Set {s[0]}..."
        )
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
        print(
            f"{pc.GREEN}fullmosaic.py  {pc.END}"
            f": Converting the mosaic to mag per square arcsec for {filter}-Band Set {s[0]}..."
        )
        psa = 2.5*n.log10((platescale[int(s[0])-1]*60)**2) # platescale adjustment
        stm1 = arcpy.sa.Raster(gridsetp + os.sep + 'skytopoc')
        stm2 = stm1 / exptime[0]
        stm3 = arcpy.sa.Log10(stm2)
        stm4 = 2.5 * stm3
        skytopomags = zeropoint[int(s[0])-1] + psa - stm4

        # Save mags mosaic to disk
        print(
            f"{pc.GREEN}fullmosaic.py  {pc.END}"
            f": Creating layer files for full-resolution mosaic for {filter}-Band Set {s[0]}..."
        )
        skytopomags.save(gridsetp + os.sep + 'skytopomags')
        layerName = dnight+'_%s_fullres%s'%(s[0],f[filter])
        layerfile = filepath.griddata+dnight+'/skytopomags%s%s.lyrx' %(f[filter],s[0])
        symbologyLayer = filepath.rasters+'magnitudes.lyrx'
        arcpy.management.MakeRasterLayer(gridsetp+'skytopomags', layerName)
        arcpy.management.ApplySymbologyFromLayer(layerName, symbologyLayer)
        arcpy.management.SaveToLayerFile(layerName, layerfile, "RELATIVE")

        # Final status update
        print(
            f"{pc.GREEN}fullmosaic.py  {pc.END}"
            f": {filter}-Band Set {s[0]} fullres mosaic COMPLETE"
        )

    
if __name__ == "__main__":
    mosaic('FCNA160803', ['1st',],'V')
    pass
