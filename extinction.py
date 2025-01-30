#-----------------------------------------------------------------------------#
#extinction.py
#
#NPS Night Skies Program
#
#Last updated: 2016/12/02
#
#This script finds the best-fit extinction coefficient and the instrumental
#zeropoint by:
#	(1) identifying the standards stars in the images 
#	(2) measuring their background-subtracted flux in [DN/s]
#	(3) computing their elevations to get airmass
#   (4) comparing the measured flux to their absolute magnitude
#   (5) Given the airmass and M-m for each star, fit for the extinction 
#       coefficient and the instrumental zeropoint
#
#Note: In order to use the ACP objects, the Python must be a 32-bit version.
#
#Input: 
#   (1) Calibrated images
#   (2) hipparcos_standards.txt
#   (3) plot_img number (optional) -- if given, the script will display the 
#       bestfit standard stars contours overlaid on the image data
#
#Output:
#   (1) extinction_stars_%s.txt -- list of the standard stars used for fitting
#   (2) extinction_fit_%s.png -- graphical display of the fitting result
#   (3) extinction_fit.txt -- best-fit extinction coefficient and zeropoint
#
#History:
#	Dan Duriscoe -- Created in visual basic as "extinction v4.vbs"
#	Li-Wei Hung -- Cleaned, improved, and translated to python
#   Zach Vanderbosch -- Updated to Python 3.11, replaced ACP/ASCOM with Astropy
#
#-----------------------------------------------------------------------------#

from astropy.io import fits
from astropy.time import Time
from glob import glob
from tqdm import trange
from scipy.optimize import curve_fit

import astropy.units as u
import astropy.wcs as wcs
import astropy.coordinates as coord
import matplotlib.pyplot as plt
import numpy as n
import pandas as pd
import warnings

# Ignore certain Astropy warnings
warnings.simplefilter('ignore', category=wcs.FITSFixedWarning)

# Local Source
from gaussian import Gaussian_2d
import filepath

#-----------------------------------------------------------------------------#

def plot_fit_result(data, x, y, popt_list):
    '''
    This module provides the visual presentation of the bestfit standard stars
    contours overlaid on the image data.
    '''
    # plot the image data
    plt.close()
    fig = plt.figure(1, figsize=(14, 5))
    ax1 = plt.subplot(121)
    im = ax1.imshow(data, interpolation='nearest', vmin=0, vmax=3000)
    fig.colorbar(im)
    plt.title('real image data')

    # add up all the contours of the standard stars in the best-fit image
    ax2 = plt.subplot(122, sharex=ax1, sharey=ax1)
    data_fitted = n.zeros_like(data)
    for i in range(len(popt_list)):
        data_fitted += Gaussian_2d((x, y), *popt_list[i]).reshape(1024, 1024)

    # display the contours on both subplots
    ax1.contour(x, y, data_fitted, 8, colors='w')
    ax2.contour(x, y, data_fitted, 8, colors='k')
    fig.colorbar(im)
    plt.title('best-fit contours over the standard stars')
    plt.show(block=False) 


def angular_separation(ra1, de1, ra2, de2):
    '''
    Compute great circle angular separation between a list
    of coordinates (ra1, de1) and a single reference
    coordinate (ra2, de2) using the Vincenty formula:
    
    https://en.wikipedia.org/wiki/Great-circle_distance

    Parameters
    ----------
    ra1: array
        An array of RA values [radians]
    de1: array
        An array of Dec values [radians]
    ra2: float
        The reference RA [radians]
    de2: float
        The reference Dec [radians]

    Returns
    -------
    sep: float
        Angular separation in degrees
    '''

    # Calculate portions of the Vincenty equation
    deltaRA = abs(ra1 - ra2)
    t1 = n.cos(de2) * n.sin(deltaRA)
    t2 = n.cos(de1) * n.sin(de2)
    t3 = n.sin(de1) * n.cos(de2) * n.cos(deltaRA)
    t4 = n.sin(de1) * n.sin(de2)
    t5 = n.cos(de1) * n.cos(de2) * n.cos(deltaRA)

    # Vincenty equation
    sep = n.arctan2(n.sqrt(t1**2 + (t2-t3)**2), t4+t5)
    sep = n.rad2deg(sep)

    return sep


def poly_sigfit(x,y,signum=5,niter=10):
    '''
    Function to perform iterative sigma clipping while
    fitting a linear trend to M-m (mag) versus airmass.

    Parameters:
    -----------
    x: array
        x values to be fit
    y: array
        y values to be fit
    signum: int
        Number of standard deviations above and below
        residuals top perform sigma clipping
    niter: int
        Number of rejection interations to perform

    Returns:
    --------
    param: list
        Best fit parameters [slope, y-intercept]
    cov: list
        Covariance matrix
    '''

    # Get copies of x/y arrays
    fit_x = n.copy(x)
    fit_y = n.copy(y)

    # Perform fits and sigma rejections
    nrej = 0
    for _ in range(niter):

        # Update number of data points being used in fit
        N = len(fit_x)

        # Try fitting linear trend to data
        param, cov = n.polyfit(fit_x, fit_y, 1, cov=True)

        # Calculate model values and residuals
        mod = n.polyval(param, fit_x)
        sigma = n.sqrt(sum((mod-fit_y)**2)/(N-1))
        residual = fit_y - mod

        # Re-define fitx/fity with sigma clipping
        fit_x = fit_x[
            (residual > -signum*sigma) & 
            (residual <  signum*sigma)
        ]
        fit_y = fit_y[
            (residual > -signum*sigma) & 
            (residual <  signum*sigma)
        ]

        # Get number of rejected points
        nrej += N - len(fit_x)

    # Get final indices of stars used for the fit
    N = len(fit_x)
    mod_fit = n.polyval(param, fit_x)
    mod_full = n.polyval(param, x)
    sigma = n.sqrt(sum((mod_fit-fit_y)**2)/(N-1))
    fit_indices = (
        (y-mod_full > -signum*sigma) & 
        (y-mod_full <  signum*sigma)
    )
    print(
        f'{nrej} out of {len(x)} stars clipped '
        f'with signum={signum} and niter={niter}.'
    )

    return param, cov, fit_indices


def extinction(dnight, sets, filter, plot_img=0):
    '''
    This module computes the extinction coefficient and the instrumental zero
    point. It returns the number of stars used for the fit and the location of
    the file containing the best-fit extinction coefficient, zeropoint, and
    their uncertainties. 
    '''

    zeropoint_dnight = []
    
    # #read in the standard star catalog
    # hips = n.loadtxt(filepath.standards+'hipparcos_standards.txt',dtype=object)
    # starn = hips[:,0]                                             #star names
    # ras, decs, v_mag, bv = n.array(hips[:,1:],dtype=n.float64).T  #star properties
    # Mag = {'V':v_mag, 'B':v_mag+bv}                    # absolute mag in V and B

    # # Convert the RA from hours to degress
    # ras = ras * 360/24

    #read inamd parse the standard star catalog
    hips = pd.read_csv(filepath.standards+'hipparcos_gaia_standards.csv')
    starn = hips['hip'].values.astype(str) # Hipparcos IDs
    ras = hips['ra'].values          # Right Ascension coords [deg]
    decs = hips['de'].values         # Declination coords [deg]
    v_mag = hips['vmag'].values      # Hipparcos V-band magnitudes [mag]
    bv = hips['b_v'].values          # Hipparcos B-V color [mag]
    Mag = {'V':v_mag, 'B':v_mag+bv}  # V and B magnitudes
    
    #define image xy coordinates
    x = n.arange(0, 1024)
    y = n.arange(0, 1024)
    x, y = n.meshgrid(x, y)
    
    #parameters specific to datasets with different filters
    k = {'V':'/', 'B':'/B/'}
    
    #loop through all the sets in that night
    for s in sets:
        calsetp = filepath.calibdata+dnight+'/S_0%s%s' %(s[0],k[filter])
        bestfit = []
        xscale = []
        yscale = []

        #read in the header to set the site object's parameter
        H = fits.getheader(calsetp+'ib001.fit',ext=0)
        site = coord.EarthLocation.from_geodetic(
            lon = H['LONGITUD']*u.deg,
            lat = H['LATITUDE']*u.deg,
            height = H['ELEVATIO']*u.m
        )
        exp = H['exptime'] #[s]
                
        # loop through each file in the set
        print(f'Processing images for Set {s[0]}...')
        images = sorted(glob(calsetp+'ib???.fit'))
        for imnum in trange(len(images)):

            # Get header and create WCS object
            fn = images[imnum]
            H = fits.getheader(fn,ext=0)
            W = wcs.WCS(H)

            # Get image observation time
            obstime = Time(H['JD'], format='jd', scale='utc')
            
            #proceed only if the plate (what plate?) is solved
            try:
                if H['PLTSOLVD']: pass
                else: continue
            except KeyError:
                continue

            # Convert RA/Dec to XY pix coords for stars covering the image
            seps = angular_separation(
                n.deg2rad(ras),
                n.deg2rad(decs),
                n.deg2rad(H['CRVAL1']),
                n.deg2rad(H['CRVAL2'])
            )
            w0 = n.where(seps < 17) # Image half diagonal = 16.97 deg
            hip_radec = [[r,d] for r,d in zip(ras[w0],decs[w0])]
            hip_xypix = W.all_world2pix(hip_radec, 0)

            # Skip the image with no nearby standard stars
            if isinstance(hip_xypix, list): 
                continue
            
            # Down select to stars strictly within the image 
            w1 = n.where(
                (hip_xypix[:,0] > 0) &
                (hip_xypix[:,0] < H['NAXIS1']) &
                (hip_xypix[:,1] > 0) &
                (hip_xypix[:,1] < H['NAXIS2'])
            )
            w1 = w0[0][w1]

            # PREVIOUS DISTANCE CONSTRAINTS
            # img_dec = abs(decs-H['CRVAL2']) < 12
            # img_ra = abs(ras-H['CRVAL1']) < (12/n.cos(n.deg2rad(decs)))
            # w1 = n.where(img_dec & img_ra)[0]   # stars 
            
            # Skip the image w/o standard stars
            if len(w1)==0:  
                continue
            
            # Get the XY pixel coordinates of the given RA/Dec locations
            radec = [[r,d] for r,d in zip(ras[w1],decs[w1])]
            xypix = W.all_world2pix(radec, 0)
            px1 = xypix[:,0]
            py1 = xypix[:,1]
            
            # Find the standard stars within 490 pixels of the image center
            # w2 = n.where(n.sqrt((px1-512)**2 + (py1-512)**2) < 490)

            # Find standard stars > buffer value from image edge
            buffer = 25  # edge buffer in pixels
            w2 = n.where(
                (px1 > 0 + buffer) &
                (px1 < 1024 - buffer) &
                (py1 > 0 + buffer) &
                (py1 < 1024 - buffer)
            )
            w3 = w1[w2]
            px, py = px1[w2], py1[w2]
            hip, ra, dec, M = starn[w3], ras[w3], decs[w3], Mag[filter][w3] 


            # Load in image data
            data = fits.getdata(fn,ext=0)
            popt_plot_list = []
            
            #fit 2D Gaussians to standard stars in the image
            for i in range(len(px)):
                #set the aperture radii
                r = ((x-px[i])**2+(y-py[i])**2)**0.5
                w = n.where(r<3)              #source aperture radius = 3 pix
                b = n.where((r>4) & (r<8))    #background aperture 4-8 pix ring
                
                #subtract background from the fitted data
                bg = n.median(data[b])
                f = data[w].ravel() - bg

                # Skip over objects near to saturation limit
                if max(f + bg) > 60000:
                    continue
                
                #fit
                guess = (px[i], py[i], 0.6, 50000)  #(x,y,std,brightness)
                with warnings.catch_warnings():
                    warnings.filterwarnings('error')
                    try:
                        popt = curve_fit(Gaussian_2d, (x[w],y[w]), f, p0=guess)[0]
                    except:
                        continue

                # Calculate elevation of the star
                star = coord.SkyCoord(
                    ra=ra[i]*u.deg, dec=dec[i]*u.deg, frame='icrs'
                )
                StarTopo = star.transform_to(
                    coord.AltAz(obstime=obstime, location=site)
                )
                elev = StarTopo.alt.deg

                # Skip stars with elevation near zero
                if elev < 5.0:
                    continue
                
                #set the acceptance threshold and record the measurement
                delta_position = n.sum(((popt-guess)**2)[0:2])   #position diff
                signal = popt[3]/bg       #brightness over the background level
                sigma = popt[2]           #sigma of the gaussian

                if sigma<2 and signal>25 and delta_position<1:
                    t = [fn[-7:],hip[i],M[i],elev]
                    t.extend(popt[:3])
                    t.append(popt[3]/H['exptime'])
                    bestfit.append(t)
                    popt_plot_list.append(popt)
            
            #plot the image overlaid with the bestfit contours
            if int(fn[-7:-4]) == plot_img:
                plot_fit_result(data, x, y, popt_plot_list)
                
            #reading in the solved plate scale is x and y image plane
            xscale.append(abs(H['CDELT1']))
            yscale.append(abs(H['CDELT2']))
            
        
        # Fit for the zeropoint and extinction coefficient
        stars = n.array((bestfit),dtype=object)
        M = n.float64(stars[:,2])            #V_mag, absolute
        elev = n.float64(stars[:,3])         #elevation[deg]
        flux = n.float64(stars[:,7])         #flux, background subtracted [DN]
        airmass = 1/n.sin(n.deg2rad(elev))   #airmass
        m = -2.5*n.log10(flux)               #v_mag, apparent
        
        # Perform fit with sigma-clipping
        param, cov, clipped_index = poly_sigfit(
            airmass, M-m, signum=5, niter=10
        )
        c, z = param                          # bestfit coefficient and zeropoint
        c_err, z_err = n.sqrt(cov.diagonal()) # uncertainties
        Nfit = sum(clipped_index)             # Number of sources used in fit

        # Save the list of stars used for calculating the zeropoint
        fmt = ['%7s','%8s','%7.2f','%9.2f','%7.1f','%6.1f','%5.2f','%7.f']
        H = 'File    Star   Magnitude Elevation   X      Y   sigma flux[DN/s]'
        fileout = filepath.calibdata+dnight+'/extinction_stars_%s_%s.txt'\
                  %(filter,s[0])
        n.savetxt(fileout,stars[clipped_index,:],fmt=fmt,header=H)
        
        sx = n.mean(xscale) * 60             #x plate scale ['/pix]
        sy = n.mean(yscale) * 60             #y plate scale ['/pix]
        sa = n.mean(xscale+yscale) * 60      #average plate scale ['/pix]
        
        fit_entry = [int(s[0]), Nfit, z, z_err, c, c_err, sx, sy, sa, exp]
        zeropoint_dnight.append(fit_entry)
                
        #plot the zeropoint and extinction coefficient fitting result
        a = n.arange(0,max(airmass),0.2)
        fig = plt.figure('zeropoint', figsize=(8,5))
        ax = fig.add_subplot(111)
        ax.plot(
            airmass[clipped_index], M[clipped_index]-m[clipped_index], 'o', 
            label=f'Hipparcos standard stars (N={Nfit})'
        )
        ax.plot(
            a, c*a+z, '-', lw=2,
            label='Best fit: %.3fx+%.3f' %(c,z)
        )
        ax.errorbar(
            0, z, z_err, fmt='o',
            label='zeropoint: %.3f+-%.3f'%(z,z_err)
        )
        ax.set_axisbelow(True)
        ax.grid(ls=':',lw=0.5,c='silver')
        ax.legend(loc=0, numpoints=1)
        ax.set_xlabel('Airmass',fontsize=14)
        ax.set_ylabel('M-m',fontsize=14)
        ax.set_title('Zeropoint and Extinction Coefficient',fontsize=14)
        imgout = filepath.calibdata+dnight+'/extinction_fit_%s_%s.png' \
                 %(filter,s[0])
        plt.savefig(imgout,dpi=200,bbox_inches='tight')
        plt.close('zeropoint')
    
    #save the bestfit zeropoint and extinction coefficient     
    fileout = filepath.calibdata+dnight+'/extinction_fit_%s.txt' %filter
    fmt = ['%4i', '%9i', '%13.3f', '%10.3f', '%12.3f', '%11.3f', '%11.3f', 
           '%7.3f', '%11.3f', '%13.1f']
    H1 = "set num_star_used zeropoint zeropoint_err extinction extinction_err "
    H2 = "x_scale y_scale avg_scale['/pix], exptime[s]"
    n.savetxt(fileout,n.array((zeropoint_dnight)),fmt=fmt,header=H1+H2)
    
    return len(stars), fileout

if __name__ == "__main__":
    pass
    extinction('SCBL170819', ['1st','2nd'], 'B')



