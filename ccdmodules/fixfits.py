#-----------------------------------------------------------------------------#
#fixfits.py
#
#NPS Night Skies Program
#
#Last updated: 2025/08/06
#
#This script provides a command line interface to fix FITS header
#values when necessary, for example if a longitude or latitude coordinate
#was entered incorrectly during data acquisition.
#
#Note: 
#
#Usage:
#
# python fixfits.py KEY VALUE [--directory, --filename]
#
# Required Arguments
# ------------------
#   KEY   : FITS header keyword to modify
#   VALUE : New value for the given FITS keyword
#
# Optional Arguments
# ------------------
#   --directory (-d)  : Directory containing FITS files 
#                       (Default = Current directory)
#   --filename  (-f)  : FITS filename including path 
#                       (Default = Find all FITS files in --directory)
#
# Example Usage to modify all FITS files in a directory:   
#
#       python fixfits.py LONGITUD -104.51142 -d G:\CCD\Data\calibdata\ROMO241004\S_01\
#
# Example Usage to modify a single FITS file:   
#
#       python fixfits.py LONGITUD -104.51142 -f G:\CCD\Data\calibdata\ROMO241004\S_01\ib001.fit
#
#Input:
#   None
#
#Output:
#   None
#
#History:
#	Zach Vanderbosch -- Created script
#
#-----------------------------------------------------------------------------#

import os
import sys
import argparse

from glob import glob
from astropy.io import fits

#-----------------------------------------------------------------------------#


def main():
    '''
    Main program for editing FITS header keyword values
    '''

    # Command line arguments
    parser = argparse.ArgumentParser()

    # Required arguments
    parser.add_argument('key', type=str, 
                        help="FITS header keyword to modify")
    parser.add_argument('value', type=str,
                        help="New value for the given FITS keyword")
    
    # Optional arguments with default values
    parser.add_argument('-d','--directory', type=str, default=".", 
                        help='Directory containing FITS files (Default = Current directory)')
    parser.add_argument('-f','--filename', type=str, default=None, 
                        help='FITS filename including path (Default = Find all FITS files in --directory)')
    args = parser.parse_args()

    # Extract arguments
    fitsDir = args.directory
    fitsName = args.filename
    fitsKey = args.key
    fitsValue = args.value

    # Get FITS files to operate on
    if fitsName is None:
        searchPath = os.path.join(fitsDir,"*.fit")
        fitsFiles = sorted(glob(searchPath))
        if len(fitsFiles) == 0:
            print(f"Did not find any FITS (.fit) files in directory {fitsDir}")
            sys.exit(1)
    else:
        if not os.path.isfile(fitsName):
            print(f"FITS file {fitsName} not found! Make sure path and filename are correct.")
            sys.exit(1)
        fitsFiles = [fitsName]

    # Make sure a FITS header keyword has been supplied
    if fitsKey is None:
        print("ERROR: A FITS header keyword must be supplied")
        print("Type 'python fixfits.py -h' for usage info")
        sys.exit(1)

    
    # Operate on each FITS file
    for f in fitsFiles:

        with fits.open(f, mode='update') as hdul:
            
            # Get the header
            hdr = hdul[0].header

            # Make sure supplied keyword exists in header
            if not fitsKey in list(hdr.keys()):
                print(f"WARNING: Keyword {fitsKey} does not exist in FITS file {f}")
                continue

            # Grab existing keyword value and match data types
            oldFitsValue = hdr[fitsKey]
            if isinstance(oldFitsValue, float):
                newFitsValue = float(fitsValue)
            elif isinstance(oldFitsValue, str):
                newFitsValue = str(fitsValue)
            elif isinstance(oldFitsValue, int):
                newFitsValue = int(fitsValue)

            # Update the value
            hdr[fitsKey] = newFitsValue
            hdul.flush()

        print(f"{f}  UPDATED")
            


# Run main during script execution
if __name__ == '__main__':
    main()
