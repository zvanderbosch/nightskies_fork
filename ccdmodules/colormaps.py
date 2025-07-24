#-----------------------------------------------------------------------------#
#colormaps.py
#
#NPS Night Skies Program
#
#Last updated: 2025/07/24
#
#This script generates a customized color map for displaying NPS night sky 
#panoramic data. This color map replicated the existing ArcGIS color map by 
#using the same RGB values. The range of the color map should be set from 14 to 
#24 magnitudes. 
#
#Input: 
#   (1) colormap_magnitudeslyr.txt
#           The RGB values from magnitudes.lyr
#
#Output:
#   (1) registered pyplot color map 'NPS_mag'
#
#History:
#	Li-Wei Hung -- Created
#
#-----------------------------------------------------------------------------#
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as n

# Local Source
import filepath  

#-----------------------------------------------------------------------------#

#RGB values originally from magnitudes.lyr
colormap_file = filepath.rasters+'colormap_magnitudeslyr.txt'
mag_start, mag_end, R, G, B = n.loadtxt(colormap_file).T    #RGB in 0-255 scale 

#RGB values in 0-1 scale
R /= 255.
G /= 255.
B /= 255.

#generate color map RGB list 
colors = []
for i in range(len(mag_start)):
    colors.append((R[i],G[i],B[i]))

#register the color map
cmap_name = 'NPS_mag'
nps_cm = mpl.colors.LinearSegmentedColormap.from_list(cmap_name, colors)
mpl.cm.register_cmap(cmap=nps_cm)

if __name__ == '__main__':
    plt.rcParams['image.cmap'] = 'NPS_mag'
    img = n.arange(10000).reshape(100,100)
    plt.imshow(img, interpolation='nearest')#, cmap=nps_mag
    cbar = plt.colorbar()
    plt.show(block=False)