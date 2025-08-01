#-----------------------------------------------------------------------------#
#progressbars.py
#
#NPS Night Skies Program
#
#Last updated: 2025/07/24
#
#This script makes an empty figure template for showing the progress of the data 
#processing stages. This template will then be returned to the main program to 
#allow the individual progress cell to be filled.  
#
#
#Note: This script was created with the inspiration from Dan Duriscoe's
#      progress_bars.xlsx.
#
#Input: 
#   None
#
#Output:
#   (1) Pop-out: A plot showing the progress of the data processing stage
#
#History:
#	Li-Wei Hung -- Created
#   Zach Vanderbosch -- Added second progress bar for process_metrics.py
#
#-----------------------------------------------------------------------------#
import matplotlib.pyplot as plt
import numpy as n
import pdb

from matplotlib.colors import PowerNorm
from numpy.random import rand

import filepath

#-----------------------------------------------------------------------------#

def bar_images(dnight, nset):
    '''
    Function that defines a progress bar figure for process_images.py
    
    Parameters:
    -----------
    dnight: array
        Array with names of data nights being processed
    nset: list
        List containing number of data sets for each data night

    Returns:
    --------
    (fig, ax): tuple
        The matplotlib figure and axes objects
    '''

    #process names
    pname = [
        'Number of sets',
        'Reduction',
        'Pointing',
        'Pointing Error',
        'Zeropoint',
        'Median Filter',
        'Coordinates',
        'Fullres Mosaic',
        'Galactic Mosaic',
        'Zodiacal Mosaic',
        'Median Mosaic'
    ] 
    
    #figure setting
    vmerge = 5  #number of row merged vertically in the first row
    hmerge = 3  #number of columns mergerd in the first column
    
    if len(dnight) < 5: nrow = 5+vmerge+1
    else: nrow = len(dnight)+vmerge+1
    
    fig = plt.figure('Progress Bar',figsize=((len(pname)+hmerge+2)/1.5,nrow*0.6))
    plt.clf()
    ax = fig.add_subplot(111)
    ax.set_facecolor('k')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    
    #make the grid
    x = n.arange(len(pname)+hmerge+2)
    y = n.arange(nrow)
                
    #Plot the gray grid lines
    [plt.plot(i*n.ones_like(y),y,'0.4',lw=1.5) for i in x[3:]] #vertical
    [plt.plot(x,i*n.ones_like(x),'0.4',lw=1.5) for i in y[vmerge:]] #horizontal
    
    #Plot "Number of sets"
    ax.text(3.5, vmerge-0.2, pname[0], color='w', horizontalalignment='center',
            verticalalignment='bottom', size='large', rotation='vertical')
            
    #Plot the process names        
    for i in range(1,len(pname)):
        ax.text(3.5+i, vmerge-0.2, pname[i], color='c', horizontalalignment='center',
                verticalalignment='bottom', size='large', weight='bold', 
                rotation='vertical')
        
    #Plot the data set names and number of sets:
    for i in range(len(dnight)):
        ax.text(1.5, vmerge+0.5+i, dnight[i], color='w', horizontalalignment='center',
                verticalalignment='center', size='large')
        ax.text(3.5, vmerge+0.5+i, nset[i], color='w', horizontalalignment='center',
                verticalalignment='center', size='large')
    
    
    #Color bar
    tmin = 0
    tmax = 10
    Bar = ax.pcolor([-1,0],[-2,-1,0],n.array(([tmin],[tmax])),cmap='Greens')
    C = plt.colorbar(Bar, orientation='horizontal',aspect=60,pad=0.05)
    C.set_label('Processing time (min)',size='medium')
    
    #more plot setting
    plt.title('', fontsize='x-large')
    plt.xlim(0,x[-1])
    plt.ylim(y[-1],0)
    
    plt.tight_layout()
    plt.draw()
    plt.ion()
    plt.pause(0.05)
    plt.show()
    
    return(fig, ax)


def bar_metrics(dnight, nset):
    '''
    Function that defines a progress bar figure for process_metrics.py
    
    Parameters:
    -----------
    dnight: array
        Array with names of data nights being processed
    nset: list
        List containing number of data sets for each data night

    Returns:
    --------
    (fig, ax): tuple
        The matplotlib figure and axes objects
    '''
    
    #process names
    pname = [
        'Number of sets',
        'All-Sources Stats',
        'Skyglow Stats',
        'Number of Stars Visible',
        'ALR Model',
        'Albedo Model',
        'Places 21K',
        'Sky Quality Indicators',
        'Make Graphics',
        'Make Summary Tables'
    ] 
    
    
    #figure setting
    vmerge = 5  #number of row merged vertically in the first row
    hmerge = 3  #number of columns merged in the first column
    
    if len(dnight) < 5: nrow = 5+vmerge+1
    else: nrow = len(dnight)+vmerge+1
    
    fig = plt.figure('Progress Bar',figsize=((len(pname)+hmerge+2)/1.5,nrow*0.6))
    plt.clf()
    ax = fig.add_subplot(111)
    ax.set_facecolor('k')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    
    #make the grid
    x = n.arange(len(pname)+hmerge+2)
    y = n.arange(nrow)
                
    #Plot the gray grid lines
    [plt.plot(i*n.ones_like(y),y,'0.4',lw=1.5) for i in x[3:]] #vertical
    [plt.plot(x,i*n.ones_like(x),'0.4',lw=1.5) for i in y[vmerge:]] #horizontal
    
    #Plot "Number of sets"
    ax.text(3.5, vmerge-0.2, pname[0], color='w', horizontalalignment='center',
            verticalalignment='bottom', size='large', rotation='vertical')
            
    #Plot the process names        
    for i in range(1,len(pname)):
        ax.text(3.5+i, vmerge-0.2, pname[i], color='c', horizontalalignment='center',
                verticalalignment='bottom', size='large', weight='bold', 
                rotation='vertical')
        
    #Plot the data set names and number of sets:
    for i in range(len(dnight)):
        ax.text(1.5, vmerge+0.5+i, dnight[i], color='w', horizontalalignment='center',
                verticalalignment='center', size='large')
        ax.text(3.5, vmerge+0.5+i, nset[i], color='w', horizontalalignment='center',
                verticalalignment='center', size='large')
    
    
    #Color bar
    tmin = 0
    tmax = 10
    Bar = ax.pcolor([-1,0],[-2,-1,0],n.array(([tmin],[tmax])),cmap='Greens')
    C = plt.colorbar(Bar, orientation='horizontal',aspect=60,pad=0.05)
    C.set_label('Processing time (min)',size='medium')
    
    #more plot setting
    plt.title('', fontsize='x-large')
    plt.xlim(0,x[-1])
    plt.ylim(y[-1],0)
    
    plt.tight_layout()
    plt.draw()
    plt.ion()
    plt.pause(0.05)
    plt.show()
    
    return(fig, ax)
    

if __name__ == "__main__":
    #dnight = []
    #d = raw_input('Please enter the dataset name: ')
    #while d != '1': 
    #    dnight.append(d)
    #    d = raw_input('Enter another dataset name. If finished, enter 1: ')
    dnight = ['FCNA160803',]
    nsets = ['1',]
    barfig, barax = bar_images(dnight, nsets)
    Z = n.loadtxt(filepath.calibdata+'FCNA160803/processtime.txt')
    M = n.ma.masked_invalid(Z)
    barax.pcolor(M,cmap='Greens',vmin=0,vmax=10)
    
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            t = Z[i,j]
            if n.isnan(t): continue
            if t < 1: texty = '%.1f' %t
            else: texty = '%d' %round(t)
            if t < 6: c = 'k'
            else: c = 'w'
            barax.text(j+0.5, i+0.5, texty, color=c, 
               horizontalalignment='center',
               verticalalignment='center', 
               size='medium')
               


