#-----------------------------------------------------------------------------#
#gaussian.py
#
#NPS Night Skies Program
#
#Last updated: 2025/07/24
#
#This script demonstrates how users can fit a 2D Gaussian to a star in an image
#that might contain more than one stars. This is a standalone script. 
#
#Input: 
#   (1) As defined in the main program 
#
#Output:
#   (1) images showing the results
#
#History:
#	Li-Wei Hung -- Created
#
#-----------------------------------------------------------------------------#
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
import numpy as n

#-----------------------------------------------------------------------------#

def Gaussian_2d(xy, x0, y0, sigma, V):
    '''
    This module returns the (x,y) value of the 2D gaussian function with the 
    given parameters. V is the volume under the curve. 

    Parameters:
    -----------
    xy: tuple
        Tuple (x,y) of x and y image pixel coordinates
    x0: float
        X-centroid of Gaussian function
    y0: float
        Y-centroid of Gaussian function
    sigma: float
        standard deviation (symmetrical in x and y directions)
    V: float
        Amplitude of Gaussian function, also the volume under the curve

    Returns:
    --------
    g1D: array
        The evaluated 2D gaussian with values flattened into a 1D-array
    '''
    x,y = xy
    g = V/(2*n.pi*sigma**2)*n.exp(-((x-x0)**2+(y-y0)**2)/(2*sigma**2))
    g1D = g.ravel()
    return g1D


if __name__ ==  "__main__": 
    '''
    An Example of fitting a 2D Gaussian to a star in an image
    '''
    # Create x and y indices
    x = n.linspace(0, 200, 201)
    y = n.linspace(0, 200, 201)
    x, y = n.meshgrid(x, y)
    
    #create stars in the data image
    star1 = Gaussian_2d((x, y), 90, 120, 5, 80000)
    star2 = Gaussian_2d((x, y), 130, 60, 5, 30000)
    star3 = Gaussian_2d((x, y), 150, 190, 5, 50000)
    
    data = star1+star2+star3
    
    # plot Gaussian_2d data generated above
    fig = plt.figure(1, figsize=(14, 5))
    ax1 = plt.subplot(121)
    im1 = ax1.imshow(data.reshape(201, 201),origin='bottom',interpolation='nearest')
    fig.colorbar(im1)
    plt.title('input data')
        
    # add some noise to the data before fitting
    data_noisy = data + 10*n.random.normal(size=data.shape)
    
    # setting the fitting aperture radius to 10 pix
    initial_guess = (85, 115, 30, 50000)
    r = ((x-initial_guess[0])**2+(y-initial_guess[1])**2)**0.5
    w = n.where(r<10)
    u = n.where(r.ravel()<10)
    
    # fit
    popt = curve_fit(Gaussian_2d,(x[w],y[w]),data_noisy[u],p0=initial_guess)[0]
    data_fitted = Gaussian_2d((x, y), *popt)

    # plot the results
    ax2 = plt.subplot(122, sharex=ax1, sharey=ax1)
    im2 = ax2.imshow(data_noisy.reshape(201, 201), origin='bottom')
    fig.colorbar(im2)    
    ax2.contour(x, y, data_fitted.reshape(201, 201), colors='w')
    plt.title('noisy data with best-fit contours')
    plt.show(block=False)


