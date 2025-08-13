#-----------------------------------------------------------------------------#
#medianfilter.py
#
#NPS Night Skies Program
#
#Last updated: 2025/07/24
#
#This script uses multiprocessing to apply median filter to each image. The 
#filter is a circle with 1 degree diameter. This filter size was selected to 
#ensure most (or all) point sources are effectively filtered out. Here, the
#Python Imaging Library (PIL/Pillow) is used to convert the FITS images to 
#TIFF images that are compatible with ArcGIS Pro to make mosaics. 
#
#Input: 
#   (1) ib###.fit
#           Calibrated image data and headers
#           (filepath.calibdata/DATANIGHT/S_0#)
#
#Output:
#   (1) median_ib###.tif
#           Median filtered images in TIFF format
#           (filepath.calibdata/DATANIGHT/S_0#/tiff)
#
#History:
#	Dan Duriscoe -- Created in java script as "med1.js"; used PixInsight
#	Li-Wei Hung -- Rewrote in python; replaced PixInsight by python functions
#   Zach Vanderbosch -- Py2 -> Py3 updates and replaced MaximDL with Astropy/PIL
#
#-----------------------------------------------------------------------------#

from astropy.io import fits
from glob import glob, iglob
from multiprocessing import Pool
from scipy.ndimage.filters import median_filter
from PIL import Image

import PIL
import itertools
import numpy as n

# Local Source
import filepath
import printcolors as pc

# Print status prefix
scriptName = 'medianfilter.py'
PREFIX = f'{pc.GREEN}{scriptName:19s}{pc.END}: '


#------------------------------------------------------------------------------#
#-------------------            Define Functions            -------------------#
#------------------------------------------------------------------------------#

def FilterImage(arg):
    '''
    Apply the median filter using the given mask size and save
    the images in .tif format.

    Parameters:
    -----------
    arg: tuple
        Tuple containing name of FITS file to be filtered, and
        a numpy array defining the median filter mask.
    '''

    # Parse input
    fn, mask = arg

    # Get FITS header/data
    with fits.open(fn) as hdul:
        fits_data = median_filter(hdul[0].data, footprint=mask)
        fits_hdr = hdul[0].header

    # Save as TIFF file. TIFF Tag IDs found here:
    # https://www.itu.int/itudoc/itu-t/com16/tiff-fx/docs/tiff6.pdf
    outTiff = 'tiff/median_%s.tif'%fn[-9:-4] #output file
    tiff_data = fits_data.astype(n.uint16)
    software_info = f'pillow (PIL) version {PIL.__version__}'
    tiff_info = {
        270: fits_hdr['OBJECT'],   # Description
        305: software_info,        # Software
        259: 1,                    # Compression (1 = None)
        282: 1058,                 # X Resolution (set to match MaximDL TIFFs)
        283: 1058,                 # Y Resolution (set to match MaximDL TIFFs)
        296: 2,                    # Resolution unit (2 = dpi)
        271: fits_hdr['INSTRUME']  # Camera Maker
    }
    tiff_output = Image.fromarray(tiff_data, mode="I;16")
    tiff_output.save(
        fn[:-9]+outTiff,
        tiffinfo = tiff_info
    )
    

#------------------------------------------------------------------------------#
#-------------------              Main Program              -------------------#
#------------------------------------------------------------------------------#

def filter(dnight, sets, filter):
    '''
    This module creats a mask and calls the FilterImage module to apply median 
    filter to the calibrated images through multiprocessing.

    Parameters:
    -----------
    dnight: str
        Name of data night to process
    sets: list
        List of data sets to process
    filter: str
        Name of photometric filter
    '''
        
    #filter paths
    F = {'V':'', 'B':'B/'}
    
    #set the mask radius to be ~0.5 degree
    # calsetp = filepath.calibdata+dnight+'/S_0%s/%s' %(sets[0][0],F[filter])
    # for fn in iglob(calsetp+'ib???.fit'):
    #     H = fits.open(fn)[0].header
    #     if 'CDELT1' in H.keys(): 
    #         plate_scale = abs(H['CDELT1']) #[deg/pix] X-axis plate scale
    #         r = n.floor(1./plate_scale/2)  #[pix] radius of the filter mask
    #         break
    
    #generate the mask
    # X, Y = n.meshgrid(n.arange(2*r+1), n.arange(2*r+1))
    # R = n.sqrt((X-r)**2+(Y-r)**2)
    # mask = n.zeros_like(R)
    # mask[n.where(R<=r)] = 1

    # Define the mask (25x25 pixels, 12.5 pixel radius or about 0.32 degrees)
    # This is the exact same mask used in the original pipeline.
    mask = n.array([
        [0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0],
        [0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0],
        [0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0],
        [0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0],
        [0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0],
        [0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0],
        [0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0],
        [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
        [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
        [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
        [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
        [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
        [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
        [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
        [0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0],
        [0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0],
        [0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0],
        [0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0],
        [0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0],
        [0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0],
        [0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0],
        [0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0]
    ])


    # Iterate through each set
    for s in sets:
        
        # Set up arguments for multiprocessed filtering
        calsetp = filepath.calibdata+dnight+'/S_0%s/%s' %(s[0],F[filter])
        images = glob(calsetp+'ib???.fit')
        arg = zip(images, itertools.repeat(mask))
        Nimages = len(images)

        # Status update
        print(f'{PREFIX}Processing images for {filter}-band Set {s[0]}...')

        # Begin filtering with multiprocessing
        threads = 5
        count = 0
        with Pool(processes=threads) as pool:
            result = pool.imap_unordered(FilterImage, arg)
            for _ in result:
                count += 1
                # Status update
                if count % 5 == 0:
                    print(
                        f'{PREFIX}{filter}-band Set {s[0]}, '
                        f'{count}/{Nimages} images complete'
                    )
        
    # Status update
    print(f'{PREFIX}{filter}-band all Sets COMPLETE')
    

    
if __name__ == "__main__":
    #import time
    #t1 = time.time()
    #filter('FCNA160803', ['1st',], 'V')
    #t2 = time.time()
    #print 'Total time: %.1f min' %((t2-t1)/60)
    pass
