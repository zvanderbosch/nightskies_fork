#-----------------------------------------------------------------------------#
#register.py
#
#NPS Night Skies Program
#
#Last updated: 2025/07/24
#
#This script will register the pointing of the image to the sky coordinates 
#using the standard stars captured in the image. Specifically, it uses the 
#PinPoiot.Plate object (through ACP) to compare the position of the stars in the
#images to the published Tycho2 catalog position. If the plate object cannot 
#solve the entire image, the input image will be cropped so the plate object 
#will only try to solve the central 200x200 pix. If the images still can't be 
#solved, it will be skipped. The solved X and Y coordinates are stored in the 
#header under 'CRVAL1' and 'CRVAL2'.
#
#Note: In order to use the ACP objects, the Python must be a 32-bit version. 
#
#Input: 
#   (1) ib###.fit
#           Calibrated images
#           (filepath.calibdata/DATANIGHT/S_0#)
#
#Output:
#   (1) ib###c.fits
#           Cropped images from the lowest altitude row used
#           for plate solving to avoid below-horizon areas.
#   (2) ib###s.fits
#           If initial solve fails, images are cropped to central
#           200x200 pixel regions and retried.
#   (3) ib###.fit
#           Original calibrated images with updated header
#
#History:
#	Dan Duriscoe -- Created in visual basic as "solve_images_v4b.vbs"
#	Li-Wei Hung -- Cleaned and translated to python
#   Zach Vanderbosch -- Updated Py2 --> Py3, removed ACP/ASCOM dependencies
#
#-----------------------------------------------------------------------------#

from glob import glob
from os.path import isfile
from astropy.io import fits
from astropy.time import Time
from multiprocessing import Pool

import os
import time
import json
import subprocess
import astropy.coordinates as coord
import astropy.units as u

# Local Source
import filepath    
import printcolors as pc

# Define print staus prefix
scriptName = 'register.py'
PREFIX = f'{pc.GREEN}{scriptName:19s}{pc.END}: '


#------------------------------------------------------------------------------#
#-------------------            Define Functions            -------------------#
#------------------------------------------------------------------------------#

def update_fits(fn, message):
    """
    Update FITS file headers with WCS solution from
    Astrometry.net for both the original image and
    the image used for solving if different from the
    original image (e.g. the masked or cropped image).

    Parameters
    ----------
    fn: str
        FITS file sent to astrometry.net
    message: list
        [solve_status (str), fn_orig (str)]
        2-element list containing the solve status
        and the name of the original FITS file that
        will be updated. The solve status can be
        'normal', 'cropped', or 'failed'.
    """

    # WCS Solution FITS header keys to save into image headers
    wcs_keys = [
        'WCSAXES', 'CTYPE1', 'CTYPE2','EQUINOX','LONPOLE','LATPOLE',
        'CRVAL1','CRVAL2','CRPIX1','CRPIX2','CUNIT1','CUNIT2',
        'CD1_1','CD1_2','CD2_1','CD2_2',
        'A_ORDER','A_0_0','A_0_1','A_0_2','A_1_0','A_1_1','A_2_0',
        'B_ORDER','B_0_0','B_0_1','B_0_2','B_1_0','B_1_1','B_2_0',
        'AP_ORDER','AP_0_0','AP_0_1','AP_0_2','AP_1_0','AP_1_1','AP_2_0',
        'BP_ORDER','BP_0_0','BP_0_1','BP_0_2','BP_1_0','BP_1_1','BP_2_0'
    ]

    # First parse the message
    solve_status = message[0]
    fn_orig = message[1]
    if solve_status == 'failed':
        with fits.open(fn_orig, mode='update') as hdul:
            hdul[0].header['PLTSOLVD'] = (False, 'Astrometric solution solved')
            hdul.flush()
        return

    # Get path to the WCS & calib files
    fn_base = fn.split("\\")[-1][:-4]
    astsetp = "{:s}/astrometry/".format(fn.split("\\")[0])
    fn_wcs = f"{astsetp}{fn_base}_wcs.fit"
    fn_calib = f"{astsetp}{fn_base}_calib.txt"

    # Load the WCS and original FITS headers
    with fits.open(fn_wcs) as hdul:
        wcs_hdr = hdul[0].header
    with fits.open(fn_orig) as hdul:
        orig_hdr = hdul[0].header

    # Create coordinate object using CRVAL values
    imgCoord = coord.SkyCoord(
        ra=wcs_hdr['CRVAL1']*u.deg,
        dec=wcs_hdr['CRVAL2']*u.deg,
        frame='icrs'
    )

    # Get pixel scale from the calibration file
    with open(fn_calib) as js:
        calib = json.load(js)
    pixscale = calib['pixscale'] / 3600. # [deg/pix] platescale

    # If the solved image is different from the original 
    # image, save the WCS solution there first
    if fn != fn_orig:
        with fits.open(fn, uint=False, mode='update') as hdul:
            H = hdul[0].header
            H['PLTSOLVD'] = (True, 'Astrometric solution solved')
            for key in wcs_keys:
                if key not in list(wcs_hdr.keys()):
                    continue
                H[key] = (wcs_hdr[key], wcs_hdr.comments[key])

            # Update the RA and DEC values in the header
            H['RA'] = imgCoord.ra.to_string(unit='hour',sep=' ',precision=2)
            H['DEC'] = imgCoord.dec.to_string(unit='deg',sep=' ',precision=1)
            
            # Add cdelt params using pixel scale value
            H.set('CDELT1', pixscale, '[deg/pixel] X-axis plate scale', before='CUNIT1')
            H.set('CDELT2', pixscale, '[deg/pixel] Y-axis plate scale', before='CUNIT1')

            # Add history
            if 'HISTORY' not in H:
                H['HISTORY'] = 'WCS created using the Astrometry.net suite'
                H['HISTORY'] = wcs_hdr['HISTORY'][1]
                H['HISTORY'] = wcs_hdr['HISTORY'][2]

            # Flush changes to file
            hdul.flush()

    # Now update the WCS header values if a cropped image was used for solving
    if solve_status == 'cropped':
        wcs_hdr['CRPIX1'] = wcs_hdr['CRPIX1'] + orig_hdr['NAXIS1']/2 - 100
        wcs_hdr['CRPIX2'] = wcs_hdr['CRPIX2'] + orig_hdr['NAXIS2']/2 - 100

    # Update original FITS file's header
    with fits.open(fn_orig, uint=False, mode='update') as hdul:
        H = hdul[0].header
        H['PLTSOLVD'] = (True, 'Astrometric solution solved')
        for key in wcs_keys:
            if key not in list(wcs_hdr.keys()):
                continue
            H[key] = (wcs_hdr[key], wcs_hdr.comments[key])

        # Update the RA and DEC values in the header
        H['RA'] = imgCoord.ra.to_string(unit='hour',sep=' ',precision=2)
        H['DEC'] = imgCoord.dec.to_string(unit='deg',sep=' ',precision=1)
        
        # Add cdelt/crota params
        H.set('CDELT1', pixscale, '[deg/pixel] X-axis plate scale', before='CUNIT1')
        H.set('CDELT2', pixscale, '[deg/pixel] Y-axis plate scale', before='CUNIT1')

        # Add history
        if 'HISTORY' not in H:
            H['HISTORY'] = 'WCS created using the Astrometry.net suite'
            H['HISTORY'] = wcs_hdr['HISTORY'][1]
            H['HISTORY'] = wcs_hdr['HISTORY'][2]

        # Flush changes to file
        hdul.flush()


def solve(fn):
    '''
    Function that sends an image to Astrometry.net for
    plate solving using the client.py script

    Parameters:
    -----------
    fn: str
        FITS filename

    Returns:
    --------
    message: list
        [solve_status (str), fn_orig (str)]
        A 2-element list containing the solve status
        and the name of the associated FITS file. The 
        solve status can be normal, cropped, or failed.
    '''
    
    fn_orig = fn
    astsetp = "%s/astrometry/"%(fn.split("\\")[0])
    m = int(fn_orig[-7:-4])

    # Get header and image coordinates
    fhdr = fits.getheader(fn,ext=0)

    # Use observation time, site location, and Alt/Az to predict RA/Dec
    obstime = Time(fhdr['JD'], format='jd', scale='utc')
    site = coord.EarthLocation.from_geodetic(
        lon = fhdr['LONGITUD']*u.deg,
        lat = fhdr['LATITUDE']*u.deg,
        height = fhdr['ELEVATIO']*u.m
    )
    imgTopoCoord = coord.SkyCoord(
        az = fhdr['AZ']*u.deg,
        alt = fhdr['ALT']*u.deg,
        obstime = obstime,
        location = site,
        frame='altaz'
    )
    fc = imgTopoCoord.transform_to(coord.ICRS())

    # Masking the area near the horizon in image 0-15
    if m < 16: 
        with fits.open(fn,uint=False) as hdul:
            f = hdul[0]
            f.data[630:] = 0.
            fn = fn[:-4]+'c.fit'
            f.writeto(fn, overwrite=True)

    #solve for astrometry; see "python client.py -help"
    fn_base = fn.split("\\")[-1][:-4]
    cmd = [
        'python', './ccdmodules/client.py', 
        '--apikey', f'{filepath.apikey}',
        '--upload', f'{fn}',
        '--parity', '1',
        '--scale-units', 'arcsecperpix',
        '--scale-est', '96.0',
        '--scale-err', '10.0', # percent
        '--ra', f'{fc.ra.deg:.6f}',
        '--dec', f'{fc.dec.deg:.6f}',
        '--radius', '12.0',
        '--corr', f'{astsetp}{fn_base}_corr.fit',
        '--calibrate', f'{astsetp}{fn_base}_calib.txt',
        '--wcs', f'{astsetp}{fn_base}_wcs.fit',
        '--solve-time', '30.0',
        '--crpix-center'
    ]
    try: 
        response = subprocess.run(cmd, timeout=None)
        response.check_returncode()
        message = ['normal', fn_orig] # files that have been solved normally
    except Exception as e:
        # trying to just solve the cropped (200x200 pix) image
        with fits.open(fn,uint=False) as hdul:
            f = hdul[0]
            l = int(len(f.data)/2)
            f.data = f.data[l-100:l+100,l-100:l+100]
            fn = fn[:-4]+'s.fit'
            f.writeto(fn, overwrite=True)

        #solve for astrometry; see "python client.py -help"
        fn_base = fn.split("\\")[-1][:-4]
        cmd = [
            'python', './ccdmodules/client.py', 
            '--apikey', f'{filepath.apikey}',
            '--upload', f'{fn}',
            '--parity', '1',
            '--scale-units', 'arcsecperpix',
            '--scale-est', '96.0',
            '--scale-err', '10.0', # percent
            '--ra', f'{fc.ra.deg:.6f}',
            '--dec', f'{fc.dec.deg:.6f}',
            '--radius', '12.0',
            '--corr', f'{astsetp}{fn_base}_corr.fit',
            '--calibrate', f'{astsetp}{fn_base}_calib.txt',
            '--wcs', f'{astsetp}{fn_base}_wcs.fit',
            '--solve-time', '30.0',
            '--crpix-center'
        ]
        try:
            response = subprocess.run(cmd, timeout=None)
            response.check_returncode()
            message = ['cropped',fn_orig] #files that have been cropped & solved
        except Exception as e:
            message = ['failed', fn_orig] #files that haven been failed to solve

    # Update the original FITS headers
    update_fits(fn, message)
        
    return message


def check_existing(fitsFile):
    '''
    Function that uses existing astrometric solutions
    when available to update FITS headers
    '''

    # Get astrometry path
    astsetp = "%s/astrometry/"%(fitsFile.split("\\")[0])

    # Copy original filename
    fn = fitsFile

    # Get image number and update first-row filenames
    m = int(fn[-7:-4])
    if m < 16:
        fn = fn[:-4]+'c.fit'

    # Look for astrometry files
    baseName = os.path.basename(fn).split(".")[0]
    calibFile = f"{astsetp}{baseName}_calib.txt"
    corrFile = f"{astsetp}{baseName}_corr.fit"
    wcsFile = f"{astsetp}{baseName}_wcs.fit"

    if isfile(calibFile) and isfile(corrFile) and isfile(wcsFile):
        # Create cropped first-row FITS file if it doesn't exist
        if not isfile(fn):
            with fits.open(fitsFile,uint=False) as hdul:
                f = hdul[0]
                f.data[630:] = 0.
                f.writeto(fn, overwrite=True)
        
        # Update FITS headers
        message = ['normal', fitsFile]
        update_fits(fn, message)
    
    # Check for cropped image solutions
    else:
        # Update the file name
        fns = fn[:-4]+'s.fit'

        # Look for astrometry files
        baseName = os.path.basename(fns).split(".")[0]
        calibFile = f"{astsetp}{baseName}_calib.txt"
        corrFile = f"{astsetp}{baseName}_corr.fit"
        wcsFile = f"{astsetp}{baseName}_wcs.fit"

        if isfile(calibFile) and isfile(corrFile) and isfile(wcsFile):

            # Create cropped first-row FITS file if it doesn't exist
            if m < 16:
                if not isfile(fn):
                    with fits.open(fitsFile,uint=False) as hdul:
                        f = hdul[0]
                        f.data[630:] = 0.
                        f.writeto(fn, overwrite=True)

            # Create cropped 200x200 pix FITS file(s) if it doesn't exist
            if not isfile(fns):
                with fits.open(fn,uint=False) as hdul:
                    f = hdul[0]
                    l = int(len(f.data)/2)
                    f.data = f.data[l-100:l+100,l-100:l+100]
                    f.writeto(fns, overwrite=True)

            # Update FITS headers
            message = ['cropped', fitsFile]
            update_fits(fns, message)
        
        # Set message for files without existing astrometric solutions
        else:
            message = ['missing', fitsFile]

    return message


#------------------------------------------------------------------------------#
#-------------------              Main Program              -------------------#
#------------------------------------------------------------------------------#

def matchstars(dnight, sets, filter, use_existing):
    '''
    Function that iterates over all calibrated images,
    sends them to Astrometry.net to be solved, and then
    updates the FITS headers with solve status and WCS
    keywords.

    Parameters:
    -----------
    dnight: str
        Name of data night to process
    sets: list
        List of data sets to process
    filter: str
        Name of photometric filter
    use_existing: bool
        Whether to use existing astrometric solutions 
        if available (TRUE) or re-do image solving (FALSE)
        

    Returns:
    --------
    (cropped_fn, failed_fn): tuple
        Tuple containing two lists of FITS file names
        that indicate which images were solved after 
        cropping to central 200x200 pixel region (cropped_fn) 
        and which images failed to solve (failed_fn).
    '''
    
    # Lists to store cropped and failed images
    cropped_fn = []
    failed_fn = []
    
    # Set number of parallel processes that can be
    # sent in to Astrometry.net
    threads = 3

    #looping through all the sets in that night
    t0 = time.time()
    for s in sets:

        # Status update
        print(f'{PREFIX}Registering images in {filter}-band Set {s[0]}...')

        # Get paths to FITS images and astrometry directory
        calsetp = filepath.calibdata + dnight + '/S_0' + s[0] + '/'
        if filter == 'V':
            astsetp = f'{calsetp}astrometry/'
        elif filter == 'B':
            astsetp = f'{calsetp}B/astrometry/'

        # Create astrometry directory if needed
        if not os.path.exists(astsetp):
            os.mkdir(astsetp)
        
        # Both V and B bands
        if filter == 'V':
            files = glob(calsetp+'ib???.fit')
        elif filter == 'B':
            files = glob(calsetp+'B/ib???.fit')

        # If use_existing = True, try updating FITS files 
        # using existing astrometric solutions
        if use_existing:

            # List to store FITS files that still need to be solved
            files_to_solve = []
            for fn in files:
                
                # Update FITS if solution already exists
                status,fitsFile = check_existing(fn)

                # Update file lists
                if status == 'normal':
                    print(f"{PREFIX}Used existing solution for {os.path.basename(fn)}")
                if status == 'cropped':
                    print(f"{PREFIX}Used existing solution for {os.path.basename(fn)}")
                    cropped_fn.append(fitsFile)
                if status == 'missing':
                    files_to_solve.append(fitsFile)

            # Solve images without existing solutions
            if len(files_to_solve) > 0:
                with Pool(processes=threads) as pool:
                    result = pool.imap_unordered(solve,files_to_solve)
                    for res in result:
                        if res[0] == 'cropped':
                            cropped_fn.append(res[1])
                        elif res[0] == 'failed':
                            failed_fn.append(res[1])

        # Solve all images from scratch
        else:
            with Pool(processes=threads) as pool:
                result = pool.imap_unordered(solve,files)
                for res in result:
                    if res[0] == 'cropped':
                        cropped_fn.append(res[1])
                    elif res[0] == 'failed':
                        failed_fn.append(res[1])

        # Sort file lists
        cropped_fn = sorted(cropped_fn)
        failed_fn = sorted(failed_fn)
        
    # Final status update
    t1 = time.time()
    print(f'{PREFIX}{filter}-band Set {s[0]} COMPLETE')
    print(f'{PREFIX}Total Solving Time = {(t1-t0)/60:.2f} minutes')

    return(cropped_fn, failed_fn)
    


if __name__ == "__main__":
    #import time
    #t1 = time.time()
    #cropped, failed = matchstars('FCNA160803', ['1st',], 'B')
    #t2 = time.time()
    #print cropped
    #print failed
    #print 'Total time: %.1f min' %((t2-t1)/60)
    pass

    
    
    
    
    
    
    