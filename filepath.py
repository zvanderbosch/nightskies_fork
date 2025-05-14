#-----------------------------------------------------------------------------#
#filepath.py
#
#NPS Night Skies Program
#
#Last updated: 2025/04/01
#
#This script defines directory file paths that point towards
#various data products needed for pipeline execution.
#
#
#Input: 
#   None
#
#Output:
#   None
#
#History:
#   Zach Vanderbosch -- Created file (2024 December)
#
#-----------------------------------------------------------------------------#


# Base directory for CCD data products
base = "C:/Users/zvanderbosch/data/CCD"

# Astrometry.net API key
apikey = 'kdvqtjbqkbkbuyzb'

# Directory where filelist.txt lives
processlist = f"{base}/Data/"

# Directory for Linearity Curve calibration files
lincurve = f"{base}/Images/Linearity Curves/"

# Directory for Flat-field calibration files
flats = f"{base}/Images/Master/"

# Directory for raw images
rawdata = f"{base}/Data/fielddata/"

# Directory for calibrated images & data products
calibdata = f"{base}/Data/calibdata/"

# Directory for mosaicked images
griddata = f"{base}/Data/griddata/"

# Directory containing TIFF world files (.tfw)
tiff = f"{base}/Data/rasters/tiff_tfws/"

# Directory containing standard star catalogs
standards = f"{base}/Data/standards/"

# Directory containing galactic/zodiacal model rasters
rasters = f"{base}/Data/rasters/"

# Directory containing ArcGIS map documents
maps = f"{base}/Data/maps/"

# Directory containing saved graphics
graphics = f"{base}/Data/graphics/"

# Directory for scripts
scripts = f"{base}/Scripts/"