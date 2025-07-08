#-----------------------------------------------------------------------------#
#pointing.py
#
#NPS Night Skies Program
#
#Last updated: 2016/11/15
#
#This script calculates the actual pointed azimuth (AZ) and altitude (ALT) using
#the solved RA and Dec values from the image headers. If the images are not 
#solved, the RA, Dec, AZ, and ALT values are interpolated.
#	(1) Read in the solved RA and Dec values from the image header.
#	(2) Update the coordinates to the observed date.
#	(3) Translate to the azimuth and altitude given the LAST and the longitude.
#   (4) Write the output to file
#   (5) Insert the interpolated AZ and ALT in the output file 
#   (6) Update the headers with the interpolated RA and Dec if the images are 
#       not solved.
#
#Note: In order to use the ACP objects, the Python must be a 32-bit version. 
#
#Input: 
#   (1) Calibrated images
#
#Output:
#   (1) pointerr_%s.txt
#
#History:
#	Dan Duriscoe -- Created in visual basic as "calc_pointing_error_v4.vbs"
#	Li-Wei Hung -- Cleaned, improved, and translated to Python
#   Davyd Betchkal -- Plotted the pointing error by image number
#   Zach Vanderbosch -- (Dec 2024) Py2 --> Py3 updates and replaced all 
#                       ACP/ASCOM commands with Astropy Time/Coords.
#
#-----------------------------------------------------------------------------#

from astropy.time import Time
from astropy.io import fits
from glob import iglob
from scipy.interpolate import UnivariateSpline

import numpy as n
import astropy.units as u
import astropy.coordinates as coord
import matplotlib.pyplot as plt
import os

# Local Source
import filepath
import ccdmodules.printcolors as pc

# Print status prefix
scriptName = 'pointing.py'
PREFIX = f'{pc.GREEN}{scriptName:19s}{pc.END}: '

#-----------------------------------------------------------------------------#
    
def interp_coord(filenames, solved_outputs):
    '''
    Interpolate the True_AZ and True_ALT for images that are not solved and 
    update the RA and DEC with the interpolated values in the header.
    '''

    # Parse the inputs
    solved, _, _, True_AZ, True_ALT = solved_outputs
    fi = n.array([int(filenames[i][-7:-4]) for i in range(len(filenames))])
    
    # The number of last image in every elevation row
    w = [0,15,30,40,45]
    
    # Iterate over elevation rows
    for i in range(len(w)-1):

        # Get indices of unsolved images within this elevation row
        wf = (fi>w[i]) & (fi<=w[i+1])
        if not any(filenames[wf]): 
            continue

        if i==0: # use the second row of images in elevation for interpolation
            wi = (solved>w[i+1]) & (solved<=w[i+2])
            k = min(3, sum(wi)-1)
            A = UnivariateSpline(solved[wi]-15, True_AZ[wi], k=1)
            E = UnivariateSpline(solved[wi]-15, True_ALT[wi]-25, k=k)
        else:
            wi = (solved>w[i]) & (solved<=w[i+1])
            k = min(3, sum(wi)-1)
            A = UnivariateSpline(solved[wi], True_AZ[wi], k=k)
            E = UnivariateSpline(solved[wi], True_ALT[wi], k=k)
            
        for fn in filenames[wf]:

            #insert the interpolated Obs_AZ and Obs_ALT 
            with fits.open(fn, mode='update') as hdul:
                H = hdul[0].header
                j = int(fn[-7:-4])
                entry = [j,H['AZ'],H['ALT'],float(A(j)),float(E(j))]
                solved_outputs = n.insert(solved_outputs,j-1,entry,axis=1)

                # Set time and site of observations
                obstime = Time(H['JD'], format='jd', scale='utc')
                site = coord.EarthLocation.from_geodetic(
                    lon = H['LONGITUD']*u.deg,
                    lat = H['LATITUDE']*u.deg,
                    height = H['ELEVATIO']*u.m
                )

                # Generate topocentric coord object
                imgTopoCoord = coord.SkyCoord(
                    az = float(A(j))*u.deg,
                    alt = float(E(j))*u.deg,
                    obstime = obstime,
                    location = site,
                    frame='altaz'
                )

                # Transform to ICRS reference frame
                imgCoord = imgTopoCoord.transform_to(coord.ICRS())
        
                # Update the RA and DEC in the header with the interpolated values
                H['RA'] = imgCoord.ra.to_string(unit='hour',sep=' ',precision=2)
                H['DEC'] = imgCoord.dec.to_string(unit='deg',sep=' ',precision=1)
                hdul.flush()
            
    return solved_outputs.T

    
def pointing_err(dnight, sets):
    '''
    This module is calculating the pointing error of each image.
    '''
    
    #looping through all the sets in that night
    for s in sets:

        # Status update
        print(f'{PREFIX}Calculating pointing error for Set {s[0]}...')

        # Define directory for calibrated FITS files
        calsetp = filepath.calibdata + dnight + '/S_0' + s[0] + '/'
        
        #read in the header to set the site object's parameter
        H = fits.open(calsetp+'ib001.fit',unit=False)[0].header
        site = coord.EarthLocation.from_geodetic(
            lon = H['LONGITUD']*u.deg,
            lat = H['LATITUDE']*u.deg,
            height = H['ELEVATIO']*u.m
        )
        
        #calculate the temperture-pressure correction for refraction
        temp = (H['AMTEMP_F']-32)/1.8 + 273                 #temperature [K]
        pres = (1-(0.0065*H['ELEVATIO']/288.15))**5.3       #pressure [atm]
        tpco = pres*(283/temp)                              #correction
        
        #refraction at 7.5 altitude
        refraction = tpco*(1/n.tan(n.deg2rad(7.5+7.31/11.9)))/60
        
        #just for V band
        solved, notsolved = [],[]
        True_AZ, True_ALT = [],[] 
        Input_AZ, Input_ALT = [],[]
        for fn in iglob(calsetp+'ib???.fit'):

            fns = fn[:-4]+'s'+fn[-4:]
            if os.path.exists(fns):
                H = fits.getheader(fns,ext=0)
            else:
                H = fits.getheader(fn,ext=0)
            
            #calculate the pointing error only if the plate is solved
            if 'PLTSOLVD' not in H:
                notsolved.append(fn)
                continue
            elif H['PLTSOLVD'] == False:
                notsolved.append(fn)
                continue
            solved.append(int(fn[-7:-4]))

            # Get the observation time
            obstime = Time(H['JD'], format='jd', scale='utc')

            # Set the RA-Dec coordinates of image center
            imgCoord = coord.SkyCoord(
                ra = H['CRVAL1']*u.deg,
                dec = H['CRVAL2']*u.deg,
                frame='icrs'
            )
            
            # Transform to AltAz topocentric reference frame
            imgTopoCoord = imgCoord.transform_to(
                coord.AltAz(obstime=obstime, location=site)
            )
            
            # Save input and true Alt/Az image coordinates
            Input_AZ.append(H['AZ'])
            Input_ALT.append(H['ALT'])
            True_AZ.append(imgTopoCoord.az.deg)

            # Correct altitude for atmospheric refraction on images 1-15
            if int(fn[-7:-4]) < 16: 
                True_ALT.append(imgTopoCoord.alt.deg + refraction)
            else:
                True_ALT.append(imgTopoCoord.alt.deg)
                          
        
        #interpolate the True_AZ for True_ALT for images that are not solved
        pterr = n.array([solved,Input_AZ,Input_ALT,True_AZ,True_ALT])
        pterr = interp_coord(n.array(notsolved), pterr)

        # calculate errors
        pErr = pterr.T
        pErr[3][n.where((pErr[1]==0)& (pErr[3]>180))] -= 360
        azmErr = (pErr[1] - pErr[3])*n.cos(n.deg2rad(pErr[4]))
        altErr = pErr[2] - pErr[4]
        totErr = n.sqrt(n.power(azmErr,2) + n.power(altErr,2))

        # Add error estimates to output
        pterr = n.concatenate(
            (pterr,
             n.reshape(azmErr,(pterr.shape[0],1)),
             n.reshape(altErr,(pterr.shape[0],1)),
             n.reshape(totErr,(pterr.shape[0],1))),
            axis=1
        )

        #create a pointing error plot
        errorPlot = plt.figure('errplot',figsize=(20,10))
        ax = errorPlot.add_subplot(111)
        plt.suptitle("Pointing Error by Image Number", fontsize=25, verticalalignment='top')
        ax.set_title("Data Set " + s[0], fontsize=20)
        ax.plot(pErr[0], azmErr, linestyle="-.", marker="o", markerfacecolor='None', 
            markersize=4, color = "darkorange", alpha=0.7, label="Azimuth Error")
        ax.plot(pErr[0], altErr, linestyle="--", marker="o", 
            markersize=4, color = "darkgreen", alpha=0.7, label="Altitude Error")
        ax.plot(pErr[0], totErr, linestyle="-", linewidth=2, marker="o", 
            markersize=6, color = "black", alpha=1, label="Total Error")
        ax.axhline(0, color="black", linestyle="-", alpha=0.5, zorder=-10)
        ax.set_ylim(-4, 4)
        ax.set_ylabel("Error in Degrees", fontsize=20, labelpad = 10)
        ax.set_xlabel("Image Number", fontsize=20, labelpad = 15)
        ax.set_xticks(n.arange(0, 50, 5))
        ax.legend(loc='upper left', markerscale=1.8, fontsize=18, framealpha=0.3)
        ax.minorticks_on()
        ax.tick_params(which='both', top=True, right=True, labelsize=15)
        ax.grid(ls=':', lw=0.5, c='silver')
        ax.text(0.5, -2.8, "Average Total Error:   " + '{:.3f}'.format(totErr.mean()) + u'\N{DEGREE SIGN}', fontsize=18)
        errorPlot.savefig(filepath.calibdata+dnight+'/pointerr_%s.png' %s[0])
        plt.close('errorplot')

        #saving the output file        
        outfile = filepath.calibdata+dnight+'/pointerr_%s.txt' %s[0]
        nformat = ['%4.f','%8.f','%8.1f','%8.2f','%8.2f','%8.2f','%8.2f','%8.2f']
        H = 'file Input_AZ Input_ALT Obs_AZ Obs_ALT   AZ_err  ALT_err  Tot_err' #column names
        n.savetxt(outfile,pterr,fmt=nformat,header=H)

        # Status update
        print(f'{PREFIX}Set {s[0]} COMPLETE')


if __name__ == "__main__":
    # pass
    print("Hi, running from the console.")
    pointing_err('FCNA160803', ['1st',])
    