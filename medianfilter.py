#-----------------------------------------------------------------------------#
#medianfilter.py
#
#NPS Night Skies Program
#
#Last updated: 2017/05/01
#
#This script uses multiprocessing to apply median filter to each image. The 
#filter is a circle with 1 degree diameter. This filter size was selected to 
#ensure most (or all) point sources are effectively filtered out. Here, MaxIM DL
#is needed to convert the fits images to tiff images that are compatible with 
#ArcGIS to make mosaics. 
#
#Input: 
#   (1) Calibrated image data
#   (2) Plate scale from image header
#
#Output:
#   (1) Median filtered images in median_ib###.tif format
#
#History:
#	Dan Duriscoe -- Created in java script as "med1.js"; used PixInsight
#	Li-Wei Hung -- Rewrote in python; replaced PixInsight by python functions
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

#-----------------------------------------------------------------------------#    
def FilterImage(arg):
    '''
    Apply the median filter using the given mask size and save the images in 
    .tif format.
    '''

    # Parse input
    fn, mask = arg

    # Print out progress messages
    m = int(fn[-7:-4])
    if m in range(0,50,5): 
        print(f'filtering images {m:d}/45')

    # Get FITS header/data
    # outFits = 'tiff/median_%s.fit'%fn[-9:-4] #temporary file
    with fits.open(fn) as hdul:
        fits_data = median_filter(hdul[0].data, footprint=mask)
        fits_hdr = hdul[0].header
        # hdul[0].data = fits_data
        # hdul.writeto(fn[:-9]+outFits,overwrite=True)

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
    


def filter(dnight, sets, filter):
    '''
    This module creats a mask and calls the FilterImage module to apply median 
    filter to the calibrated images through multiprocessing.
    '''
        
    #filter paths
    F = {'V':'', 'B':'B/'}
    
    #set the mask radius to be ~0.5 degree
    calsetp = filepath.calibdata+dnight+'/S_0%s/%s' %(sets[0][0],F[filter])
    for fn in iglob(calsetp+'ib???.fit'):
        H = fits.open(fn)[0].header
        if 'CDELT1' in H.keys(): 
            plate_scale = abs(H['CDELT1']) #[deg/pix] X-axis plate scale
            r = n.floor(1./plate_scale/2)  #[pix] radius of the filter mask
            break
    
    #generate the mask
    X, Y = n.meshgrid(n.arange(2*r+1), n.arange(2*r+1))
    R = n.sqrt((X-r)**2+(Y-r)**2)
    mask = n.zeros_like(R)
    mask[n.where(R<=r)] = 1
    
    #multiprocessing the filtering process
    p = Pool()
    
    #loop through all the sets in that night
    for s in sets:
        calsetp = filepath.calibdata+dnight+'/S_0%s/%s' %(s[0],F[filter])
        arg = zip(glob(calsetp+'ib???.fit'), itertools.repeat(mask))
        p.map(FilterImage, arg)
        
    p.close()
    p.join()
    

    
if __name__ == "__main__":
    #import time
    #t1 = time.time()
    #filter('FCNA160803', ['1st',], 'V')
    #t2 = time.time()
    #print 'Total time: %.1f min' %((t2-t1)/60)
    pass
