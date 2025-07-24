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
#           (filepath.flats)
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

def reducev(dnight, sets, flatname, curve):
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
    '''
    #read in the linearity curve (ADU, multiplying factor)
    xp, fp = n.loadtxt(filepath.lincurve+curve+'.txt', unpack=True, delimiter=",")
    
    #read in the flat
    flat = fits.open(filepath.flats+flatname,unit=False)[0].data
    
    #looping through all the sets in that night
    for s in sets:
        rawsetp = filepath.rawdata + dnight + '/' + s + '/'
        calsetp = filepath.calibdata + dnight + '/S_0' + s[0] + '/'
        print(f'{PREFIX}Reducing {dnight} V-band Set {s[0]}...')
        if os.path.isdir(calsetp):
            print(f'{PREFIX}Replacing old calibrated files...')
        else:
            os.makedirs(calsetp)
            os.makedirs(calsetp+'tiff/')
    
        #correct the darks for bias and linearity response; crop the biases
        dark = []
        bias = []
        for i in range(1,6):
            darkraw = fits.open(rawsetp+'dark%i.fit'%i,uint=False)[0]
            biasraw = fits.open(rawsetp+'mbias%i.fit'%i,uint=False)[0]
            biasrawd = biasraw.data
            darkp = (darkraw.data - biasrawd).clip(0) #replace negatives with 0
            darki = darkp * n.interp(darkp,xp,fp) #correct linearity response
            biascrop = biasrawd[486:536,486:536]
            dark.append(darki)
            bias.append(biasrawd)
            fits.writeto(calsetp+'thermal%i.fit'%i, darki, overwrite=True)
            fits.writeto(rawsetp+'biasc%i.fit'%i, biascrop, overwrite=True,
                         header=biasraw.header)
        
        #average combine to generate the master thermal and bias
        corthermal = n.average(dark,axis=0)
        combias = n.average(bias,axis=0)
        fits.writeto(calsetp+'corthermal.fit', corthermal, overwrite=True)
        fits.writeto(calsetp+'combias.fit', combias, overwrite=True)
        
        #measure the bias drift for each frames
        nb = len(glob(rawsetp+'biasc*.fit'))
        biasc = n.empty([nb,50,50])
        Temp = n.empty([nb,1])                    # CCD temperature [C]
        for i in range(1,nb+1):
            biasci = fits.open(rawsetp+'biasc%i.fit'%i,uint=False)[0]
            biasc[i-1] = biasci.data
            Temp[i-1] = biasci.header['CCD-TEMP']
        baseline = n.average(biasc[:6])
        biasdrift_full = n.average(biasc,axis=(1,2)) - baseline
        biasdrift = biasdrift_full[5:]
        n.savetxt(filepath.calibdata+dnight+'/biasdrift_%s.txt'%s[0],
                  biasdrift,fmt='%5.3f',header='delta_bias[ADU]')
        fig = plt.figure('bias')
        plt.plot(n.arange(len(biasc)), n.zeros(len(biasc)), 'k--')
        plt.plot(n.arange(len(biasc)), biasdrift_full, 'o')
        plt.ylim(-5,5)
        plt.title('Bias Drift Compared to the Average of the First 5 Files')
        plt.xlabel('Bias File number')
        plt.ylabel('Delta_Bias [ADU] (bias - %i)'%baseline)
        plt.savefig(filepath.calibdata+dnight+'/biasdrift_%s.png'%s[0])   
        plt.close('bias')
        
        #calibrate the science images
        file = n.hstack((rawsetp+'zenith1.fit',
                         glob(rawsetp+'ib???.fit'),
                         rawsetp+'zenith2.fit'))
            
        for i in range(len(file)):

            with fits.open(file[i],uint=False) as hdu:
                f = hdu[0]                            # science image 
                f.data -= combias+biasdrift[i]        # subtract drift-corrected bias
                f.data *= n.interp(f.data,xp,fp)      # correct for linearity response
                f.data -= corthermal                  # subtract dark
                f.data /= flat                        # divide by flat
                f.data = f.data.clip(min=1.0)         # Set minimum value to 1
                f.data = f.data.clip(max=65535.)      # Set maximum value to saturation limit
                f.data = f.data.astype(n.uint16)      # Convert to uint16 values
                f.header['IMAGETYP'] = 'CALIB_M'      # Update header
                f.writeto(calsetp+file[i][len(rawsetp):], overwrite=True)

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
                calsetp+'tiff/'+file[i][len(rawsetp):-4]+'.tif',
                tiffinfo = tiff_info
            )
        
        for f in iglob(filepath.tiff+'*.tfw'):
            shutil.copy2(f,calsetp+'tiff/')
                    

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
    flat = fits.open(filepath.flats+flatname,unit=False)[0].data

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
        
        for f in iglob(filepath.tiff+'*.tfw'):
            shutil.copy2(f,calsetp+'tiff/')
                        