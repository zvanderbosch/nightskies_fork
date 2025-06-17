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
from scipy.optimize import curve_fit
from photutils.aperture import aperture_photometry
from photutils.aperture import CircularAperture
from photutils.aperture import CircularAnnulus
from photutils.aperture import ApertureStats
from astropy.stats import SigmaClip
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Circle,Annulus
from astropy.visualization import ZScaleInterval

import os
import stat
import shutil
import astropy.units as u
import astropy.wcs as wcs
import astropy.coordinates as coord
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as n
import pandas as pd
import warnings


# Ignore certain Astropy warnings
warnings.simplefilter('ignore', category=wcs.FITSFixedWarning)

# Local Source
from gaussian import Gaussian_2d
import printcolors as pc
import filepath

# Print status prefix
scriptName = 'extinction.py'
PREFIX = f'{pc.GREEN}{scriptName:19s}{pc.END}: '

ZS = ZScaleInterval(
    n_samples=10000, contrast=0.15, max_reject=0.5, 
    min_npixels=5, krej=2.5, max_iterations=5
)
mpl.rcParams.update(
    {# Use mathtext, not LaTeX
    'text.usetex': False,
    'axes.formatter.use_mathtext': True,
    'font.family': 'STIXGeneral',
    'mathtext.fontset': 'cm',
    'axes.unicode_minus': False
    }
)

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


def linear_func_Zfree(x,slope,yint):
    '''Linear model with slope and y-intercept as free parameters'''
    return slope*x + yint


def linear_func_Zfixed(x,slope):
    '''Linear model with only slope as free parameter'''
    return slope*x


def poly_sigfit(x,y,signum=5,niter=10,fixedZ=False):
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


        # Try fitting linear trend to data and calculate model
        if fixedZ:
            pguess = [-0.2]
            param, cov = curve_fit(linear_func_Zfixed, fit_x, fit_y, p0=pguess)
            mod = linear_func_Zfixed(fit_x, *param)
        else:
            pguess = [-0.2, 14.7]
            param, cov = curve_fit(linear_func_Zfree, fit_x, fit_y, p0=pguess)
            mod = linear_func_Zfree(fit_x, *param)

        # Calculate residuals
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
    if fixedZ:
        mod_fit = linear_func_Zfixed(fit_x, *param)
        mod_full = linear_func_Zfixed(x, *param)
    else:
        mod_fit = linear_func_Zfree(fit_x, *param)
        mod_full = linear_func_Zfree(x, *param)
    sigma = n.sqrt(sum((mod_fit-fit_y)**2)/(N-1))
    fit_indices = (
        (y-mod_full > -signum*sigma) & 
        (y-mod_full <  signum*sigma)
    )

    return param, cov, fit_indices


def remove_readonly(func, path, excinfo):
    '''
    Error-catching function to handle removal of read-only folders
    '''
    os.chmod(path, stat.S_IWRITE)
    func(path)


def extinction(dnight, sets, filter, zeropoint, plot_img=0):
    '''
    This module computes the extinction coefficient and the instrumental zero
    point. It returns the number of stars used for the fit and the location of
    the file containing the best-fit extinction coefficient, zeropoint, and
    their uncertainties. 
    '''

    # List to save fit results for each dataset
    zeropoint_dnight = []

    # Convert deafult zeropoint to float
    zeropoint = float(zeropoint)
    
    #read in the standard star catalog
    # hips = n.loadtxt(filepath.standards+'hipparcos_standards.txt',dtype=object)
    # starn = hips[:,0]                                             #star names
    # ras, decs, v_mag, bv = n.array(hips[:,1:],dtype=n.float64).T  #star properties
    # Mag = {'V':v_mag, 'B':v_mag+bv}                    # absolute mag in V and B

    # # Convert the RA from hours to degress
    # ras = ras * 360/24

    #read inamd parse the standard star catalog
    # hips = pd.read_csv(filepath.standards+'hipparcos_gaia_standards.csv')
    hips = pd.read_csv(filepath.standards+'hipparcos_gaia_standards_6pixAper.csv')
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

        # Check for figure directory if image cutout plotting is enabled.
        # Plotting will usually be disabled except for diagnostic pusposes.
        makePlots = False
        if makePlots:
            figsetp = f"{calsetp}cutouts/"
            if os.path.exists(figsetp):
                shutil.rmtree(figsetp, onerror=remove_readonly)
            os.makedirs(figsetp)
                
        # loop through each file in the set
        print(f'{PREFIX}Processing images for {filter}-band Set {s[0]}...')
        images = sorted(glob(calsetp+'ib???.fit'))
        imagesSolved = 0
        for imnum in range(len(images)):

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
            imagesSolved += 1

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
            
            # Skip the image w/o standard stars
            if len(w1)==0:
                continue
            
            # Get the XY pixel coordinates of the given RA/Dec locations
            radec = [[r,d] for r,d in zip(ras[w1],decs[w1])]
            xypix = W.all_world2pix(radec, 0)
            px1 = xypix[:,0]
            py1 = xypix[:,1]
            
            # Find the standard stars within 490 pixels of the image center
            w2 = n.where(n.sqrt((px1-512)**2 + (py1-512)**2) < 490)

            # Get final selection of stars
            w3 = w1[w2]
            px, py = px1[w2], py1[w2]
            hip, ra, dec, M, bvcol = starn[w3], ras[w3], decs[w3], Mag[filter][w3] , bv[w3]

            # Load in image data
            data = fits.getdata(fn,ext=0)
            popt_plot_list = []
            
            #fit 2D Gaussians to standard stars in the image
            xposImg, yposImg = [], []
            for i in range(len(px)):

                # Skip stars with elevation below 12-degrees
                star = coord.SkyCoord(
                    ra=ra[i]*u.deg, dec=dec[i]*u.deg, frame='icrs'
                )
                StarTopo = star.transform_to(
                    coord.AltAz(obstime=obstime, location=site)
                )
                elev = StarTopo.alt.deg
                if elev < 12.0:
                    continue

                #set the aperture radii
                apRadius = 6  # Dan's Pipeline uses  6, this pipeline originally used 3
                annInner = 8  # Dan's Pipeline uses  8, this pipeline originally used 4
                annOuter = 12 # Dan's Pipeline uses 12, this pipeline originally used 8
                r = ((x-px[i])**2+(y-py[i])**2)**0.5
                w = n.where(r<apRadius)
                b = n.where(
                    (r>annInner) & 
                    (r<annOuter)
                )
                
                #subtract background from the fitted data
                bg = n.median(data[b])
                f = data[w].ravel() - bg

                # Skip over objects near to saturation limit
                if (max(f + bg) > 55000) | (max(f + bg) < 4000):
                    continue
                
                # Perform initial Gaussian fit
                guess = (px[i], py[i], 0.6, 50000)  #(x,y,std,brightness)
                with warnings.catch_warnings():
                    warnings.filterwarnings('error')
                    try:
                        popt = curve_fit(Gaussian_2d, (x[w],y[w]), f, p0=guess)[0]
                        pxr,pyr = popt[0],popt[1] # Refined XY centroid
                    except:
                        continue

                # Perform a refined Gaussian fit
                r = ((x-pxr)**2+(y-pyr)**2)**0.5
                w = n.where(r < apRadius)
                b = n.where((r > annInner) & (r < annOuter))
                bg = n.median(data[b])
                f = data[w].ravel() - bg
                guess = (pxr, pyr, popt[2], popt[3])  #(x,y,std,brightness)
                with warnings.catch_warnings():
                    warnings.filterwarnings('error')
                    try:
                        popt = curve_fit(Gaussian_2d, (x[w],y[w]), f, p0=guess)[0]
                        pxf,pyf = popt[0],popt[1] # Final XY centroid
                    except:
                        continue

                # Perform aperture photometry with photutils
                position = [(pxr,pyr)] # Using "refined" XY centroid to match PSF fitting
                aperture = CircularAperture(position, r=apRadius)
                annulus = CircularAnnulus(
                    position, r_in=annInner, r_out=annOuter
                )
                photTable = aperture_photometry(data, aperture)

                # Subtract median background in annulus from central aperture
                sigclip = SigmaClip(sigma=5.0, maxiters=10)
                bkgStats = ApertureStats(data, annulus, sigma_clip=sigclip)
                bkgSum = bkgStats.median * aperture.area
                photTable['Flux'] = (photTable['aperture_sum'] - bkgSum) / H['exptime']
                
                #set the acceptance threshold and record the measurement
                # delta_position = n.sum(((popt-guess)**2)[0:2])   #position diff
                # signal = popt[3]/bg       #brightness over the background level
                # sigma = popt[2]           #sigma of the gaussian
                # if sigma<2 and signal>25 and delta_position<2:

                # Save result as long as flux measurements are positive
                if (popt[3] > 0) & (photTable['Flux'] > 0):
                    t = [fn[-7:],hip[i],M[i],elev,bvcol[i]]
                    t.extend(popt[:3])
                    t.append(popt[3]/H['exptime'])
                    t.append(photTable['Flux'])
                    bestfit.append(t)
                    popt_plot_list.append(popt)
                    xposImg.append(popt[0])
                    yposImg.append(popt[1])

                    # Make summary figure showing image cutouts
                    if makePlots:

                        fig = plt.figure(figsize=(17,5))
                        gs = GridSpec(1,3)
                        ax = fig.add_subplot(gs[0])
                        bx = fig.add_subplot(gs[1])
                        cx = fig.add_subplot(gs[2])

                        # Get area of image to plot
                        imw = 14
                        implot = data[
                            int(py[i]-imw):int(py[i]+imw),
                            int(px[i]-imw):int(px[i]+imw)
                        ]
                        vmin,vmax = ZS.get_limits(implot)

                        # Plot image in first two panels
                        ax.imshow(data,vmin=vmin,vmax=vmax,cmap='Greys')
                        bx.imshow(data,vmin=vmin,vmax=vmax,cmap='Greys')

                        # Add initial guess and fitted centroid
                        ax.scatter(px[i],py[i],fc='c',ec='k',s=30)
                        bx.scatter(px[i],py[i],fc='c',ec='k',s=30,label='Initial Guess')
                        ax.scatter(pxr,pyr,fc='gold',ec='gold',marker='+',s=80)
                        bx.scatter(pxr,pyr,fc='gold',ec='gold',marker='+',s=80,label='Refined Centroid')
                        ax.scatter(pxf,pyf,fc='r',ec='r',marker='x',s=80)
                        bx.scatter(pxf,pyf,fc='r',ec='r',marker='x',s=80,label='Final Centroid')

                        # Add legends
                        bx.legend(loc='upper left',framealpha=0.8)

                        # Add apertures and annuli
                        rt = 10. / n.sqrt(2.) # 10 over root 2
                        circAper1 = Circle((pxr,pyr), 6.0, ec='forestgreen', fc='None', lw=2)
                        circAper2 = Circle((pxr,pyr), apRadius, ec='forestgreen', fc='None', lw=2)
                        annulus1 = Annulus((pxr,pyr), 12.0, 4.0, ec='m', fc='None', lw=1.5, alpha=0.5, ls=':')
                        annulus2 = Annulus((pxr,pyr), annOuter, annOuter-annInner, ec='m', fc='m', lw=2, alpha=0.2, ls='--')
                        aperAnnulus1 = Circle((pxr+10, pyr+ 0), 2.0, ec='m', fc='m', alpha=0.2, lw=1.5)
                        aperAnnulus2 = Circle((pxr+ 0, pyr+10), 2.0, ec='m', fc='m', alpha=0.2, lw=1.5)
                        aperAnnulus3 = Circle((pxr-10, pyr+ 0), 2.0, ec='m', fc='m', alpha=0.2, lw=1.5)
                        aperAnnulus4 = Circle((pxr+ 0, pyr-10), 2.0, ec='m', fc='m', alpha=0.2, lw=1.5)
                        aperAnnulus5 = Circle((pxr+rt, pyr+rt), 2.0, ec='m', fc='m', alpha=0.2, lw=1.5)
                        aperAnnulus6 = Circle((pxr-rt, pyr+rt), 2.0, ec='m', fc='m', alpha=0.2, lw=1.5)
                        aperAnnulus7 = Circle((pxr-rt, pyr-rt), 2.0, ec='m', fc='m', alpha=0.2, lw=1.5)
                        aperAnnulus8 = Circle((pxr+rt, pyr-rt), 2.0, ec='m', fc='m', alpha=0.2, lw=1.5)
                        ax.add_patch(circAper1)
                        ax.add_patch(annulus1)
                        ax.add_patch(aperAnnulus1)
                        ax.add_patch(aperAnnulus2)
                        ax.add_patch(aperAnnulus3)
                        ax.add_patch(aperAnnulus4)
                        ax.add_patch(aperAnnulus5)
                        ax.add_patch(aperAnnulus6)
                        ax.add_patch(aperAnnulus7)
                        ax.add_patch(aperAnnulus8)
                        bx.add_patch(circAper2)
                        bx.add_patch(annulus2)

                        # Plot Radial Profile in third panel
                        ropt = ((x-pxf)**2+(y-pyf)**2)**0.5
                        rmod = n.arange(0.0,7.025,0.025)
                        gmod = popt[3]/(2*n.pi*popt[2]**2)*n.exp(-(rmod**2)/(2*popt[2]**2))
                        cx.scatter(ropt[w],f/1e3,c='k')
                        cx.axhline(0.0,ls=':',c='k',lw=1)
                        cx.plot(rmod, gmod/1e3, c='r')

                        # Set XY limits
                        ax.set_xlim(pxr-imw, pxr+imw)
                        ax.set_ylim(pyr-imw, pyr+imw)
                        bx.set_xlim(pxr-imw, pxr+imw)
                        bx.set_ylim(pyr-imw, pyr+imw)
                        cx.set_xlim(0,apRadius)

                        # Set XY labels
                        ax.set_xlabel('X (pix)',fontsize=16)
                        bx.set_xlabel('X (pix)',fontsize=16)
                        cx.set_xlabel('Radial Distance (pix)',fontsize=16)
                        ax.set_ylabel('Y (pix)',fontsize=16)
                        bx.set_ylabel('Y (pix)',fontsize=16)
                        cx.set_ylabel('Flux (10$^3$ DN)',fontsize=16)

                        # Tick params
                        ax.tick_params(which='both',labelsize=14)
                        bx.tick_params(which='both',labelsize=14)
                        cx.tick_params(which='both',labelsize=14)

                        # Add titles
                        suptitle = f"HIP {hip[i]} in Image {fn[-9:]}, {dnight}"
                        titlea = 'Original Pipeline Apertures'
                        titleb = 'Updated Pipeline Apertures'
                        titlec = 'Updated Pipeline Radial Profile'
                        plt.suptitle(suptitle,fontsize=20,y=1.0)
                        ax.set_title(titlea, fontsize=15)
                        bx.set_title(titleb, fontsize=15)
                        cx.set_title(titlec, fontsize=15)

                        figname = f'{figsetp}cutout_{hip[i]}_aper{apRadius:.0f}.png'
                        plt.savefig(figname, dpi=200, bbox_inches='tight')
                        plt.close()
            
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
        #flux = n.float64(stars[:,8])        #PSF flux, background subtracted [DN]
        flux = n.float64(stars[:,9])         #Aperture flux, background subtracted [DN]
        m = -2.5*n.log10(flux)               #v_mag, apparent
        mdiffRaw = M - m                     #Magnitude difference (absolute - apparent)

        # Use Hardie (1962) equation for airmass
        za = 90. - elev
        secantZ = 1./n.cos(n.deg2rad(za))
        am1 = 0.0018167 * (secantZ - 1.0)
        am2 = 0.0028750 * (secantZ - 1.0)**2
        am3 = 0.0008083 * (secantZ - 1.0)**3
        airmass = secantZ - am1 - am2 - am3

        # Apply distortion correction (?? not really sure what this is about)
        xpos = stars[:,5].astype(n.float64)
        ypos = stars[:,6].astype(n.float64)
        offcenter = n.sqrt((xpos-511)**2 + (ypos-511)**2)
        mcorr = n.zeros_like(m, dtype=n.float64)
        mcorr[offcenter > 340] = 0.0003 * (offcenter[offcenter > 340] - 340)
        m += mcorr

        # Apply color correction with fixed color coefficient
        colorCoeff = 0.04
        color = stars[:,4].astype(n.float64)
        mdiff = M - m
        mdiff += colorCoeff * color
        
        # Perform fit with zeropoint and extinction as free parameters
        paramFree, covFree, clipped_index = poly_sigfit(
            airmass, mdiff, signum=5, niter=10, fixedZ=False
        )
        # Perform fit with fixed default zeropoint
        paramFixed, covFixed, _ = poly_sigfit(
            airmass, mdiff-zeropoint, signum=5, niter=10, fixedZ=True
        )
        # Perform fit with zeropoint and color-coefficient as free parameters
        cparam, ccov, _ = poly_sigfit(
            color, mdiffRaw, signum=5, niter=10, fixedZ=False
        )

        # Extract Parameters
        extFree, zpFree = paramFree                           # bestfit extinction and zeropoint
        extFixed = paramFixed[0]                              # bestfit extinction with fixed zeropoint
        colorCoeffFree = cparam[0]                            # bestfit color coefficient
        extFree_err, zpFree_err = n.sqrt(covFree.diagonal())  # uncertainties
        extFixed_err = n.sqrt(covFixed.diagonal())[0]         # uncertainties
        colorCoeffFree_err = n.sqrt(ccov.diagonal())[0]       # uncertainties
        Nfit = sum(clipped_index)                             # Number of sources used in fit
        Nrej = sum(~clipped_index)                            # Number of sources rejected in fit


        # Save the list of stars used for calculating the zeropoint
        fmt = ['%7s','%8s','%7.2f','%9.2f','%7.2f','%7.1f','%6.1f','%5.2f','%8.f','%8.f']
        H = 'File    Star   Magnitude Elevation   B-V   X      Y      sigma   PSF_flux  Aper_flux'
        fileout = filepath.calibdata+dnight+'/extinction_stars_%s_%s.txt'\
                  %(filter,s[0])
        n.savetxt(fileout,stars[clipped_index,:],fmt=fmt,header=H)
        
        # Calculate average plate scales
        sx = n.mean(xscale) * 60             #x plate scale ['/pix]
        sy = n.mean(yscale) * 60             #y plate scale ['/pix]
        sa = n.mean(xscale+yscale) * 60      #average plate scale ['/pix]
        
        # Save fit results to list
        fit_entry = [
            int(s[0]), imagesSolved, Nfit, Nrej,
            zpFree, zpFree_err, extFree, extFree_err, 
            zeropoint, extFixed, extFixed_err,
            colorCoeff, colorCoeffFree, colorCoeffFree_err,
            sx, sy, sa, exp
        ]
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
            a, extFree*a + zpFree, '-', lw=2, c='C1',
            label='Best fit (ZP-Free): %.3fx+%.3f' %(extFree,zpFree)
        )
        ax.plot(
            a, extFixed*a + zeropoint, '--', lw=2, c='k',
            label='Best fit (ZP-Fixed): %.3fx+%.2f' %(extFixed,zeropoint)
        )
        ax.errorbar(
            0, zpFree, zpFree_err, fmt='o', mfc='None', mec='C1',
            label='Free Zeropoint: %.3f+-%.3f'%(zpFree,zpFree_err)
        )
        ax.errorbar(
            0, zeropoint, 0, fmt='o', mfc='None', mec='k', ecolor='None',
            label='Fixed Zeropoint: %.2f'%(zeropoint)
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

        # Status update
        print(f'{PREFIX}{filter}-band Set {s[0]} COMPLETE')
    
    #save the bestfit zeropoint and extinction coefficient     
    fileout = filepath.calibdata+dnight+'/extinction_fit_%s.txt' %filter
    fmt = [
        '%5i'   , '%11i'  , '%14i'  , '%18i'  , # H1 columns
        '%15.3f', '%19.3f', '%16.3f', '%20.3f', # H2 columns
        '%18.2f', '%18.3f', '%22.3f',           # H3 columns
        '%20.3f', '%17.3f', '%21.3f',           # H4 columns
        '%8.3f' , '%8.3f' , '%17.3f', '%11.1f'  # H5 columns
    ]
    H1 = "set  img_solved  num_star_used  num_star_rejected  "
    H2 = "zeropoint_free  zeropoint_free_err  extinction_free  extinction_free_err  "
    H3 = "zeropoint_default  extinction_fixedZ  extinction_fixedZ_err  "
    H4 = "color_coeff_default  color_coeff_free  color_coeff_free_err  "
    H5 = "x_scale  y_scale  avg_scale['/pix]  exptime[s]"
    n.savetxt(fileout,n.array((zeropoint_dnight)),fmt=fmt,header=H1+H2+H3+H4+H5)
    
    return len(stars), fileout

if __name__ == "__main__":
    pass
    extinction('SCBL170819', ['1st','2nd'], 'B')



