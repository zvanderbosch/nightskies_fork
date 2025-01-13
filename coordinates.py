#-----------------------------------------------------------------------------#
#cooordinates.py
#
#NPS Night Skies Program
#
#Last updated: 2016/12/16
#
#This script calculates ecliptic and galactic coordinates and rotation angles of
#the images. The output will be used in producing natural sky model for a given 
#location, date, and time.
#
#Input: 
#   (1) Calibrated and solved image headers
#   (2) Obs_AZ and Obs_ALT from pointerr_%s.txt
#
#Output:
#   (1) coordinates_%s.txt with all units in degree
#
#History:
#	Dan Duriscoe -- Created in visual basic "compute image coordinates v4.vbs"
#	Li-Wei Hung -- Cleaned, improved, and translated to python
#
#-----------------------------------------------------------------------------#

from astropy import units as u
from astropy.io import fits
from astropy.time import Time
from glob import glob

import numpy as n
import astropy.coordinates as coord

# Local Source
import filepath

#-----------------------------------------------------------------------------#
def bearing_angle(lat1, lon1, lat2, lon2):
    '''
    Calculate the bearing angle of (lat2, lon2) with respect to (lat1, lon1). 
    The bearing angle ranges from 0 to 360 degrees, with zero at due north 
    and increasing clockwise. Both the inputs and outputs are in degrees. 
    '''
    lat1, lon1, lat2, lon2 = n.deg2rad([lat1, lon1, lat2, lon2])
    x = n.cos(lat1)*n.sin(lat2) - n.sin(lat1)*n.cos(lat2)*n.cos(lon1-lon2)
    y = n.sin(lon1-lon2)*n.cos(lat2)
    bearing = n.rad2deg(n.arctan(y/x))
    if x < 0: bearing += 180
    return -bearing % 360.


def wcs_position_angle(hdr):
    '''
    Calculates position angle of image using WCS CD matrix.
    Using code from Astrometry.net (see get_orientation):
    https://github.com/dstndstn/astrometry.net/blob/main/net/wcs.py

    Parameters:
    -----------
    hdr: Astropy FITS Header
        FITS header with WCS CD matrix values

    Returns:
    --------
    orient: float
        Position angle in range 0-360, measured CCW east from north
    '''

    # Find parity of image
    detCD = hdr['CD1_1']*hdr['CD2_2'] - \
            hdr['CD1_2']*hdr['CD2_1']
    if detCD >= 0:
        parity = 1.
    else:
        parity = -1.
    
    # Get positive position angle (0-360), measured CCW East from North
    T = parity * hdr['CD1_1'] + hdr['CD2_2']
    A = parity * hdr['CD2_1'] - hdr['CD1_2']
    orient = -n.rad2deg(n.arctan2(A,T)) + 180.

    return orient

    
def galactic_ecliptic_coords(dnight, sets):
    '''
    This module computes the galactic and ecliptic coordinates needed for 
    building the natural sky model. 
    '''
    
    # Define the ecliptic and galactic N-poles in RA-Dec coords
    ecliptic_pole = [66.56, 18.]            #N pole Dec [deg] and RA [hr]
    galactic_pole = [27.1283, 167.1405]     #N pole latitude and longitude [deg]


    #loop through all the sets in that night
    for s in sets:
        calsetp = filepath.calibdata+dnight+'/S_0%s/' %s[0]
        outlist = []
        
        #read in the header to set the site object's parameter
        H = fits.getheader(calsetp+'ib001.fit',ext=0)
        site = coord.EarthLocation.from_geodetic(
            lon = H['LONGITUD']*u.deg,
            lat = H['LATITUDE']*u.deg,
            height = H['ELEVATIO']*u.m
        )
        
        #read in the registered images coordinates
        fnum, Obs_AZ, Obs_ALT = n.loadtxt(
            f'{filepath.calibdata}{dnight}/pointerr_{s[0]}.txt', 
            usecols=(0,3,4)
        ).T

        # loop through each file in the set
        file_list = sorted(glob(calsetp+'ib???.fit'))
        for i,fn in enumerate(file_list):

            H = fits.getheader(fn,ext=0)
            w = n.where(fnum==int(fn[-7:-4]))

            # Get the observation time from header
            obstime = Time(H['JD'], format='jd', scale='utc')

            #------------- Calculate the galactic coordinates
            ra, dec = H['RA'], H['DEC']
            c = coord.SkyCoord(
                ra, dec, 
                unit=(u.hourangle, u.deg), 
                distance=100*u.kpc,
                equinox=obstime
            )

            # Convert image center RA/Dec to Galactic coords
            galactic_l = -c.galactic.l.wrap_at(180*u.deg).deg
            galactic_b = c.galactic.b.deg
            
            if ('PLTSOLVD' in H.keys()) and H['PLTSOLVD']:   #if plate is solved
                pa = wcs_position_angle(H)
                b = bearing_angle(c.dec.deg, -c.ra.deg, *galactic_pole)
                galactic_angle = (b + pa) % 360
            else: 
                # Generate topocentric coord object
                galPoleCoord = coord.SkyCoord(
                    ra = galactic_pole[1]*u.deg,
                    dec = galactic_pole[0]*u.deg,
                    frame='icrs'
                )
                galPoleTopo = galPoleCoord.transform_to(
                    coord.AltAz(obstime=obstime, location=site)
                )
                galactic_angle = bearing_angle(
                    Obs_ALT[w][0], 
                    Obs_AZ[w][0], 
                    galPoleTopo.alt.deg, 
                    galPoleTopo.az.deg
                )

            #------------- Calculate the ecliptic coordinates

            eclPoleCoord = coord.SkyCoord(
                ra = ecliptic_pole[1]*360/24*u.deg,
                dec = ecliptic_pole[0]*u.deg,
                frame='icrs'
            )
            eclPoleTopo = eclPoleCoord.transform_to(
                coord.AltAz(obstime=obstime, location=site)
            )
            ecliptic_angle = bearing_angle(
                Obs_ALT[w][0], 
                Obs_AZ[w][0], 
                eclPoleTopo.alt.deg, 
                eclPoleTopo.az.deg
            )

            csun = coord.get_sun(obstime)
            ecliptic_l = -(c.heliocentrictrueecliptic.lon.degree-csun.ra.degree)
            ecliptic_l = (ecliptic_l+180)%360-180
            ecliptic_b = c.heliocentrictrueecliptic.lat.degree
            
            outlist.append([int(fn[-7:-4]), galactic_angle, galactic_l, 
            galactic_b , ecliptic_angle, ecliptic_l, ecliptic_b, c.ra.deg, c.dec.deg]) #[deg]
            
        #save the coordinates
        fmt = ['%5i','%8.2f','%8.2f','%8.2f','%8.2f','%8.2f','%8.2f','%12.6f','%10.6f']
        H = 'File  Gal_ang   Gal_l    Gal_b   Ecl_ang   Ecl_l    Ecl_b   RA          Dec      [deg]'
        fileout = filepath.calibdata+dnight+'/coordinates_%s.txt'%s[0]
        n.savetxt(fileout,n.array(outlist),fmt=fmt,header=H)


if __name__ == "__main__":
    #galactic_ecliptic_coords('FCNA160803', ['1st',])
    pass
