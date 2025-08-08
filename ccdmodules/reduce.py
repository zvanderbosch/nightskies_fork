#-----------------------------------------------------------------------------#
#reduce.py
#
#NPS Night Skies Program
#
#Last updated: 2025/07/24
#
#This script performs basic image reduction on the image data collected by the 
#NPS Night Skies Program. The script corrects images for:
#	(1) Bias 
#	(2) Dark 
#	(3) Flat
#   (4) Linearity response of the detector
#
#Note: Camera temperature must be matched to linearity curve
#
#Input:
#   (1) ib###.fit, zenith#.fit, biasc#.fit, mbias#.fit, dark#.fit
#           Raw dark, bias, and science frames 
#           (filepath.fielddata/DATANIGHT)
#   (2) Master Flat FITS file
#           Master flat-field calibration file 
#           (filepath.calimages)
#   (3) Linearity curve TXT file
#           Linearity curve calibration file 
#           (filepath.lincurve)
#
#Output:
#   (1) combias.fit
#           Master bias calibration file
#           (filepath.calibdata/DATANIGHT/S_0#)
#   (2) corthermal.fit
#           Master dark calibration file
#           (filepath.calibdata/DATANIGHT/S_0#)
#   (3) ib###.fit
#           Calibrated images in FITS format
#           (filepath.calibdata/DATANIGHT/S_0#)
#   (4) ib###.tif
#           Calibrated images in TIFF format
#           (filepath.calibdata/DATANIGHT/S_0#/tiff)
#   (5) ib###.tfw, median_ib###.tif
#           TIFF World files (.tfw), copied from filepath.rasters/tiff_tfws
#           (filepath.calibdata/DATANIGHT/S_0#/tiff)
#   (5) biasdrift_<DATASET>.txt
#           Measured bias drift per image in ADU
#           (filepath.calibdata/DATANIGHT)
#   (6) biasdrift_<DATASET>.png
#           Figure showing brias drift data per image
#           (filepath.calibdata/DATANIGHT)
#
#History:
#	Dan Duriscoe -- Created in 2011 in visual basic as "calibrate images.vbs"
#	Li-Wei Hung -- Cleaned and translated to python
#   Zach Vanderbosch -- Updated Py2 -> Py3 and removed MaximDL dependencies
#
#-----------------------------------------------------------------------------#

from astropy.io import fits
from glob import glob, iglob
from PIL import Image

import PIL
import matplotlib.pyplot as plt
import numpy as n
import os
import shutil

# Local Source
import filepath
import printcolors as pc

# Define print staus prefix
scriptName = 'reduce.py'
PREFIX = f'{pc.GREEN}{scriptName:19s}{pc.END}: '

#------------------------------------------------------------------------------#
#-------------------            Define Functions            -------------------#
#------------------------------------------------------------------------------#

def reducev(dnight, sets, flatname, curve, standards):
    '''
    This module is for calibrating the V-band data.

    Parameters:
    -----------
    dnight: string
        Name of data night to process (e.g. ROMO241004)
    sets: list
        List of data sets to process
    flatname: string
        Name of master flat field FITS file
    curve: string
        Name of linearity curve TXT file
    standards: bool
        Whether to use standard dark/bias frames located
        in Images/Master (TRUE) or use the dark/bias 
        frames taken during data collection (FALSE).
    '''

    # Get filenames for standard calibration files
    curveFile = f"{filepath.lincurve}{curve}.txt"
    flatFile = f"{filepath.calimages}{flatname}"
    biasFile = f"{filepath.calimages}{curve}bias.fit"
    darksecFile = f"{filepath.calimages}{curve}thermal_1sec.fit"

    # Read in the linearity curve (ADU, multiplying factor)
    xp, fp = n.loadtxt(curveFile, unpack=True, delimiter=",")
    
    # Read in the standard flat, dark, and bias frames
    flat = fits.getdata(flatFile, ext=0, unit=False)
    biasStd = fits.getdata(biasFile, ext=0, unit=False)
    darksecStd = fits.getdata(darksecFile, ext=0, unit=False)
    
    # Loop through all sets in the given night
    for s in sets:
        
        # Status update
        print(f'{PREFIX}Reducing {dnight} V-band Set {s[0]}...')

        # Define raw and calibrated data directories
        rawsetp = f"{filepath.rawdata}{dnight}/{s}/"
        calsetp = f"{filepath.calibdata}{dnight}/S_0{s[0]}/"

        # Get filenames for science images
        scienceFiles = n.hstack(
            (
                rawsetp+'zenith1.fit',
                sorted(glob(rawsetp+'ib???.fit')),
                rawsetp+'zenith2.fit'
            )
        )
        Nf = len(scienceFiles)
        
        # Check that calibdata directory exists
        if os.path.isdir(calsetp):
            print(f'{PREFIX}Replacing old calibrated files...')
        else:
            os.makedirs(calsetp)
            os.makedirs(calsetp+'tiff/')

        # Use standard dark/bias calibration frames
        if standards:

            # Get exposure time
            with fits.open(f"{rawsetp}zenith1.fit") as hdul:
                hdr = hdul[0].header
                exptime = hdr['EXPTIME']

            # Scale 1-second dark exposure to exposure time
            darkStd = (darksecStd * exptime).clip(0)

            # Perform linearity correction
            darkStdCorr = darkStd * n.interp(darkStd,xp,fp)

            # Define the master thermal, bias, and cropped bias images
            corthermal = darkStdCorr
            combias = biasStd
            biascrop = combias[486:536,486:536]

        # Use dark/bias frames from data collection night
        else:

            # Load dark/bias frames and correct the darks 
            # for bias and linearity response
            dark = []
            bias = []
            for i in range(5):
                darkraw = fits.getdata(f"{rawsetp}dark{i+1}.fit", ext=0, uint=False)
                biasraw = fits.getdata(f"{rawsetp}mbias{i+1}.fit", ext=0, uint=False)
                darkp = (darkraw - biasraw).clip(0) #replace negatives with 0
                darki = darkp * n.interp(darkp,xp,fp) #correct linearity response
                dark.append(darki)
                bias.append(biasraw)
            
            # Average combine to generate master thermal, bias, and cropped bias images
            corthermal = n.average(dark,axis=0)
            combias = n.average(bias,axis=0)
            biascrop = combias[486:536,486:536]


        # Save master thermal, bias, and cropped bias to files
        fits.writeto(f"{calsetp}corthermal.fit", corthermal, overwrite=True)
        fits.writeto(f"{calsetp}combias.fit", combias, overwrite=True)
        fits.writeto(f"{calsetp}combiasc.fit", biascrop, overwrite=True)
        
        # Measure the bias drift for each frame
        baseline = n.average(biascrop)
        biasc = n.empty([Nf,50,50])
        for i in range(Nf):
            biasc[i] = fits.getdata(
                f"{rawsetp}biasc{i+6}.fit", ext=0, uint=False
            )
        biasdrift = n.average(biasc,axis=(1,2)) - baseline

        # Save bias drift data to text file
        n.savetxt(
            f"{filepath.calibdata}{dnight}/biasdrift_{s[0]}.txt",
            biasdrift, fmt='%5.3f', header='delta_bias[ADU]'
        )
        
        # Create bias drift figure
        fig = plt.figure('bias')
        ax = fig.add_subplot(111)
        ax.plot(n.arange(Nf), n.zeros(len(biasc)), 'k--')
        ax.plot(n.arange(Nf), biasdrift, 'o')
        ax.set_ylim(-5,5)
        ax.set_title('Bias Drift Compared to the Average of the First 5 Files')
        ax.set_xlabel('Bias File number')
        ax.set_ylabel('Delta_Bias [ADU] (bias - %i)'%baseline)
        plt.savefig(f"{filepath.calibdata}{dnight}/biasdrift_{s[0]}.png")   
        plt.close('bias')
        
        # Perform image bias/dark/flat calibration  
        for i,fn in enumerate(scienceFiles):

            # Get base filename
            fnBase = os.path.basename(fn)

            with fits.open(fn,uint=False) as hdu:
                f = hdu[0]                            # science image 
                f.data -= combias+biasdrift[i]        # subtract drift-corrected bias
                f.data *= n.interp(f.data,xp,fp)      # correct for linearity response
                f.data -= corthermal                  # subtract dark
                f.data /= flat                        # divide by flat
                f.data = f.data.clip(min=1.0)         # Set minimum value to 1
                f.data = f.data.clip(max=65535.)      # Set maximum value to saturation limit
                f.data = f.data.astype(n.uint16)      # Convert to uint16 values
                f.header['IMAGETYP'] = 'CALIB_M'      # Update header
                f.writeto(f"{calsetp}{fnBase}", overwrite=True)

            # Save as TIFF file. TIFF Tag IDs found here:
            # https://www.itu.int/itudoc/itu-t/com16/tiff-fx/docs/tiff6.pdf
            tiff_data = f.data.astype(n.uint16)
            software_info = f'pillow (PIL) version {PIL.__version__}'
            tiff_info = {
                270: f.header['OBJECT'],   # Description
                305: software_info,        # Software
                259: 1,                    # Compression (1 = None)
                282: 1058,                 # X Resolution (set to match MaximDL TIFFs)
                283: 1058,                 # Y Resolution (set to match MaximDL TIFFs)
                296: 2,                    # Resolution unit (2 = dpi)
                271: f.header['INSTRUME']  # Camera Maker
            }
            tiff_output = Image.fromarray(tiff_data, mode="I;16")
            tiff_output.save(
                f"{calsetp}tiff/{fnBase[:-4]}.tif",
                tiffinfo = tiff_info
            )
        
        # Copy TIFF world files (.tfw) to calibdata tiff directory
        for f in iglob(f"{filepath.tiff}*.tfw"):
            shutil.copy2(f,f"{calsetp}tiff/")
                    

def reduceb(dnight, sets, flatname, curve):
    '''
    This module is for calibrating the B-band data. Some of the computation is 
    dependent from the output from the reducev module.

    Parameters:
    -----------
    dnight: string
        Name of data night to process (e.g. ROMO241004)
    sets: list
        List of data sets to process
    flatname: string
        Name of master flat field FITS file
    curve: string
        Name of linearity curve TXT file
    '''
    #read in the linearity curve (ADU, multiplying factor)
    xp, fp = n.loadtxt(filepath.lincurve+curve+'.txt', unpack=True, delimiter=",")
    
    #read in the flat
    flat = fits.open(filepath.calimages+flatname,unit=False)[0].data

    #looping through all the sets in that night
    for s in sets:
        rawsetp = filepath.rawdata + dnight + '/' + s + '/'
        calsetp = filepath.calibdata + dnight + '/S_0' + s[0] + '/B/'
        print(f'{PREFIX}Reducing {dnight} B-band Set {s[0]}...')
        if os.path.isdir(calsetp):
            print(f'{PREFIX}Replacing old calibrated files...')
        else:
            os.makedirs(calsetp)
            os.makedirs(calsetp+'tiff/')
    
        #read in the thermal, bias, and bias drift from the V band directory
        corthermal = fits.open(calsetp[:-2]+'corthermal.fit')[0].data * 1.5
        combias = fits.open(calsetp[:-2]+'combias.fit')[0].data
        biasdrift = n.loadtxt(filepath.calibdata+dnight+'/biasdrift_%s.txt'%s[0], unpack=True)[1:-1]
        
        #calibrate the science images
        file = glob(rawsetp+'ib???b.fit')
        
        for i in range(len(file)):
            f = fits.open(file[i],uint=False)[0]  # science image
            f.data -= combias+biasdrift[i]        # subtract drift-corrected bias
            f.data *= n.interp(f.data,xp,fp)      # correct for linearity response
            f.data -= corthermal                  # subtract dark
            f.data /= flat                        # divide by flat
            f.data = f.data.clip(min=1.0)         # Set minimum value to 1
            f.data = f.data.astype(n.uint16)      # Convert to uint16 values
            f.header['IMAGETYP'] = 'CALIB_M'
            f.writeto(calsetp+file[i][len(rawsetp):-5]+'.fit', overwrite=True)

            # Save as TIFF file. TIFF Tag IDs found here:
            # https://www.itu.int/itudoc/itu-t/com16/tiff-fx/docs/tiff6.pdf
            tiff_data = f.data.astype(n.uint16)
            software_info = f'pillow (PIL) version {PIL.__version__}'
            tiff_info = {
                270: f.header['OBJECT'],   # Description
                305: software_info,        # Software
                259: 1,                    # Compression (1 = None)
                282: 1058,                 # X Resolution (set to match MaximDL TIFFs)
                283: 1058,                 # Y Resolution (set to match MaximDL TIFFs)
                296: 2,                    # Resolution unit (2 = dpi)
                271: f.header['INSTRUME']  # Camera Maker
            }
            tiff_output = Image.fromarray(tiff_data, mode="I;16")
            tiff_output.save(
                calsetp+'tiff/'+file[i][len(rawsetp):-5]+'.tif',
                tiffinfo = tiff_info
            )
        
        # Copy TIFF world files (.tfw) to calibdata tiff directory
        for f in iglob(filepath.tiff+'*.tfw'):
            shutil.copy2(f,calsetp+'tiff/')
                        