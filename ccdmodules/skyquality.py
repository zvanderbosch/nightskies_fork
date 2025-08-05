#-----------------------------------------------------------------------------#
#skyquality.py
#
#NPS Night Skies Program
#
#Last updated: 2025/07/25
#
#This script computes sky quality index (SQI) and synthetic
#sky quality meter (SQM) values, along with a few photometric
#indicators using simple sky brightness aperture photometry
#measurements.
#
#Note: 
#
#Input:
#   (1) extinction_fit_<FILTER>.xlsx
#           Best-fit extinction parameters
#           (filepath.calibdata/DATANIGHT)
#   (2) ib###.fit, zenith#.fit
#           Calibrated images
#           (filepath.calibdata/DATANIGHT/S_0#)
#   (3) sqitbl.df, sqitbl80.dbf, sqitbl70.dbf
#           Sky quality index (SQI) histogram tables from skyglow module
#           (filepath.griddata/DATANIGHT/S_0#)
#   (4) calibreport.xlsx
#           Calibration report spreadsheet
#           (filepath.calibdata/DATANIGHT)
#   (5) skybright_positions.xlsx
#           Spreadsheet containing pixel locations per image
#           where photometric measurements will be made
#           (filepath.spreadsheets)
#
#Output:
#   None
#
#
#History:
#	Zach Vanderbosch -- Created script (translated from sqiv4.vbs)
#
#-----------------------------------------------------------------------------#

from astropy.io import fits
from dbfread import DBF
from photutils.aperture import aperture_photometry
from photutils.aperture import CircularAperture
from photutils.aperture import ApertureStats
from astropy.stats import SigmaClip
from astropy.time import Time

import os
import numpy as n
import pandas as pd
import astropy.units as u
import astropy.coordinates as coord

# Local Source
import filepath
import printcolors as pc

# Define print staus prefix
scriptName = 'skyquality.py'
PREFIX = f'{pc.GREEN}{scriptName:19s}{pc.END}: '

#------------------------------------------------------------------------------#
#-------------------            Define Functions            -------------------#
#------------------------------------------------------------------------------#

def get_site_info(imageFile):
    '''
    Function to grab observing site and time info

    Parameters:
    -----------
    imageFile: str
        Path to FITS file used for site information

    Returns:
    --------
    obsTime: astropy.time.Time
        Time object representation observation time in UTC
    obsSite: EarthLocation 
        Astropy EarthLocation object for observing site
    '''

    # Load image header
    H = fits.getheader(imageFile)

    # Set up site and time astropy objects
    obsTime = Time(
        H['DATE-OBS'], format='isot', scale='utc'
    )
    obsSite = coord.EarthLocation.from_geodetic(
        lon = H['LONGITUD']*u.deg,
        lat = H['LATITUDE']*u.deg,
        height = H['ELEVATIO']*u.m
    )

    return obsTime, obsSite


def measure_skybrightness(imgPath):
    '''
    Function to perform sky brightness aperture
    photometry for SQM calculations.

    Parameters:
    -----------
    imgPath: str
        Path to directory where calibrated FITS files exist

    Returns:
    --------
    allPhot: DataFrame
        Dataframe containing photometric measurements
    '''

    # Load in skybright positions spreadsheet
    allPositions = pd.read_excel(
        f"{filepath.spreadsheets}skybright_positions.xlsx"
    )

    # Iterate over each mosaic image
    allPhot = []
    for i in range(45):

        # Get image data
        imgFile = f"{imgPath}ib{i+1:03d}.fit"
        if not os.path.isfile(imgFile):
            continue
        with fits.open(imgFile) as hdul:
            imgData = hdul[0].data

        # Get pixel positions for given image
        imgPositions = allPositions.loc[allPositions.image == i+1]
        xpix = imgPositions.PixelX.values
        ypix = imgPositions.PixelY.values
        xyPositions = [(x,y) for x,y in zip(xpix,ypix)]

        # Perform aperture photometry on image
        apRadius = 20.
        apertures = CircularAperture(xyPositions, r=apRadius)
        photTable = aperture_photometry(imgData, apertures)
        photTable = photTable.to_pandas()

        # Calculate aperture median
        sigclip = SigmaClip(sigma=5.0, maxiters=10)
        aperStats = ApertureStats(imgData, apertures, sigma_clip=sigclip)
        imgPositions.insert(2,'ADU',aperStats.median)

        # Save photometry to list
        allPhot.append(imgPositions)

    # Concatenate photometry into single dataframe
    allPhot = pd.concat(allPhot,ignore_index=True)
    
    return allPhot


def measure_zenith_images(imgPath):
    '''
    Function to perform sky brightness aperture
    photometry for just the central portion of zenith images.

    Parameters:
    -----------
    imgPath: str
        Path to directory where calibrated FITS files exist

    Returns:
    --------
    allPhot: DataFrame
        Dataframe containing photometric measurements
    '''

    # Measure sky brightness at center of zenith images
    allPhot = []
    for i in range(2):
    
        # Get image data
        imgFile = f"{imgPath}zenith{i+1}.fit"
        if not os.path.isfile(imgFile):
            continue
        with fits.open(imgFile) as hdul:
            imgData = hdul[0].data
            imgHdr = hdul[0].header

        # Perform aperture photometry on image
        apRadius = 20.
        apertures = CircularAperture([(511,511)], r=apRadius)
        photTable = aperture_photometry(imgData, apertures)
        photTable = photTable.to_pandas()

        # Calculate aperture median
        sigclip = SigmaClip(sigma=5.0, maxiters=10)
        aperStats = ApertureStats(imgData, apertures, sigma_clip=sigclip)
        photTable['aperture_median'] = aperStats.median
        photTable['inst_mag'] = 2.5 * n.log10(photTable['aperture_median'] / imgHdr['EXPTIME'])
        photTable['image'] = [f'zenith{i+1}']
        allPhot.append(photTable)

    # Concatenate photometry into single dataframe
    allPhot = pd.concat(allPhot,ignore_index=True)
    
    return allPhot


def planet_brightness(planet, time, site, dEarthSun):
    '''
    Function to calculate apparent magnitude of a
    solar system planet at a given time and location.

    Parameters:
    -----------
    planet: str
        Name of planet to process
    time: astropy.time.Time
        Observation time (UTC)
    site: EarthLocation
        Astropy EarthLocation for observation site
    dEarthSun: float
        Distance between Earth and Sun at given time (AU)

    Returns:
    --------
    pmag: float
        Planet apparent magnitude (mag)
    pTopoCoord.az: float
        Planet Azimuth coordinate (deg)
    pTopoCoord.alt: float
        Planet Altitude coordinate (deg)
    '''

    # Metadata for each planet
    planetMetadata = {
        'venus': {
            'albedo': 0.00002,
            'exponent': 0.65,
            'zeropoint': 27.5
        },
        'mars': {
            'albedo': 0.0000085,
            'exponent': 0.5,
            'zeropoint': 26.7
        },
        'jupiter': {
            'albedo': 0.0002,
            'exponent': 0.65,
            'zeropoint': 27.9
        },
        'saturn': {
            'albedo': 0.000174,
            'exponent': 0.5,
            'zeropoint': 27.75
        }
    }
    albedo = planetMetadata[planet]['albedo']
    exp = planetMetadata[planet]['exponent']
    zp = planetMetadata[planet]['zeropoint']

    # Get topocentric and heliocentric coordinates
    pCoord = coord.get_body(planet, time, site)
    pTopoCoord = pCoord.transform_to(
        coord.AltAz(obstime=time, location=site)
    )
    pHelioCoord = pCoord.transform_to(
        coord.HeliocentricMeanEcliptic()
    )

    # Get distances to earth and sun in Astronomical Units
    distToEarth = pTopoCoord.distance.value
    distToSun = pHelioCoord.distance.value

    # Calculate flux factor for earth-sun-planet positioning
    t1 = distToEarth**2 + distToSun**2 - dEarthSun**2
    t2 = 2 * distToEarth * distToSun
    fluxFactor = (1 + t1/t2) / 2

    # Calculate apparent magnitude
    t1 = distToEarth * distToSun
    t2 = albedo * fluxFactor**exp
    pmag = 5.0 * n.log10(t1/t2) - zp
    
    return pmag, pTopoCoord.az.deg, pTopoCoord.alt.deg


def calculate_airmass(za):
    '''
    Function to calculate airmass using Kasten & Young (1989) equation

    Parameters:
    -----------
    za: array
        Zenith angles (deg)
    
    Returns:
    --------
    airmass: array
        Airmass values
    '''

    # Use Kasten & Young (1989) equation for airmass
    cosZ = n.cos(n.deg2rad(za))
    am1 = 0.50572 * (96.07995 - za)**(-1.6364)
    airmass = 1. / (cosZ + am1)

    return airmass


def calc_SQI(gridPath, mask):
    '''
    Function to calculate Sky Quality Index (SQI)

    Parameters:
    -----------
    gridPath: str
        Path to data night's griddata directory
    mask: str
        String indicating the mask to use. Must be horizon, ZA80, or ZA70.

    Returns:
    --------
    sqi: float
        The sky quality index value
    '''

    # Set static weighting factors
    weightingNaturalSky = n.array(
        [0,2,3.2,5,7.9,10,10,10,10,10,10,10,10,10,
         10,10,10,10,10,10,10,10,10,10,10,10,10,10],
        dtype=n.float64
    )
    weightingMilkyWay = n.array(
        [0,1.3,2.1,3.3,5.2,6.65,7.05,7.3,7.7,8.35,9.3,10,
         10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10],
        dtype=n.float64
    )
    weightingScotopicVision = n.array(
        [0,0,0,0,0,0,0,0,0,0,0.55,1.1,1.65,2.2,2.75,3.3,
         3.85,4.4,4.95,5.5,6.05,6.6,7.15,7.7,8.25,8.8,9.35,10],
        dtype=n.float64
    )
    weightingStarsVisible = n.array(
        [0,0.1,0.25,0.45,0.87,1.25,1.7,2.5,3.75,4.85,5.8,6.75,7.35,
         7.85,8.2,8.5,8.8,9,9.2,9.35,9.5,9.6,9.7,9.8,9.87,9.92,9.97,10], 
        dtype=n.float64
    )

    # Load in SQI table corresponding to mask
    if mask == 'horizon':
        sqiTableFile = f"{gridPath}sqitbl.dbf"
    elif mask == 'ZA80':
        sqiTableFile = f"{gridPath}sqitbl80.dbf"
    elif mask == 'ZA70':
        sqiTableFile = f"{gridPath}sqitbl70.dbf"
    sqiTableDBF = DBF(sqiTableFile)
    sqiTable = pd.DataFrame(iter(sqiTableDBF))

    # Get histogram values and order from darkest to brightest bins
    histValues = sqiTable['Value_0'].values[::-1]
    sumValues = sqiTable['Value_0'].sum()

    # Calculate frequency of values in each bin
    frequency = 100 * histValues / sumValues

    # Calculate component linear index values
    indexNaturalSky = 0.25 * sum(weightingNaturalSky * frequency) / 10
    indexMilkyWay = 0.25 * sum(weightingMilkyWay * frequency) / 10
    indexStarsVisible = 0.25 * sum(weightingStarsVisible * frequency) / 10

    # Calculate component log index for scotopic vision
    linearScotopicSum = sum(weightingScotopicVision * frequency)
    if linearScotopicSum <= 0:
        logScotopicSum = 0.0
    else:
        logScotopicSum = n.log10(linearScotopicSum + 1.0)
    logFactor = 100. / n.log10(1000.)
    indexScotopicVision = 0.25 * logFactor * logScotopicSum

    # Calculate final hybrid SQI
    sqi = 100 - indexNaturalSky - indexMilkyWay - indexScotopicVision - indexStarsVisible

    return sqi    


def calc_skyonly_SQM(dataNight, setNum, filterName):
    '''
    Function to calculate synthetic Sky Quality Meter (SQM)
    metric using just the measured background sky brightness.

    Parameters:
    -----------
    dataNight: str
        Name of data night
    setNum: int
        Data set number
    filterName: str
        Name of photometric filter

    Returns:
    --------
    sqm: float
        Synthethic sky-only Sky Quality Meter (SQM) value
    photometry: DataFrame
        DataFrame containing photometric measurements
    '''
    
    # Get zeropoint, extinction coeff, plate scale, & exposure time
    extFile = f"{filepath.calibdata}{dataNight}/extinction_fit_{filterName}.xlsx"
    extData = pd.read_excel(extFile)
    zeropoint = extData['zeropoint_default'].iloc[setNum-1]
    platescale = extData['avg_scale'].iloc[setNum-1]
    exptime = extData['exptime'].iloc[setNum-1]
    psa = 2.5*n.log10((platescale*60)**2) # platescale adjustment

    # Perform sky brightness measurements
    imgsetp = f"{filepath.calibdata}{dataNight}/S_{setNum:02d}/"
    photometry = measure_skybrightness(imgsetp)

    # Calculate cosine of Zenith angle and apply to ADU
    photometry['cosZ'] = n.cos(n.deg2rad(90.0 - photometry.Pany.values))
    photometry['net_arcsec'] = photometry.cosZ * 3600 * 3600 * 4
    photometry['cosADU'] = photometry.cosZ * photometry.ADU

    # Get sums for zenith angles <= 54 degrees
    arcsecSum = photometry[photometry.Pany >= 36].net_arcsec.sum()
    aduSum = photometry[photometry.Pany >= 36].cosADU.sum()

    # Calculate SQM scale factor
    sqmScaleFactor = 2.5 * n.log10(arcsecSum)

    # Convert platescale adjustment (mag) to sq. arcsec per ADU
    psaADU = (4*3600*3600) / (10**(0.4*psa))

    # Calculate synthetic SQM
    mags = zeropoint - 2.5 * n.log10(aduSum * psaADU / exptime)
    sqm = sqmScaleFactor + mags

    return sqm, photometry


def calc_synthetic_SQM(dataNight, setNum, filterName, skySQM):
    '''
    Function to calculate synthetic Sky Quality Meter (SQM)
    metric using measured background sky brightness and
    including flux contributed from visible stars and planets.

    Parameters:
    -----------
    dataNight: str
        Name of data night
    setNum: int
        Data set number
    filterName: str
        Name of photometric filter
    skySQM: float
        Synthetic sky-only SQM from calc_skyonly_SQM

    Returns:
    --------
    sqm: float
        Synthethic Sky Quality Meter (SQM) including sky and planets
    '''

    # Get zeropoint, extinction coeff, plate scale, & exposure time
    extFile = f"{filepath.calibdata}{dataNight}/extinction_fit_{filterName}.xlsx"
    extData = pd.read_excel(extFile)
    extCoeff = abs(extData['extinction_fixedZ'].iloc[setNum-1])

    # Get Site time and location for data set midpoint
    imgsetp = f"{filepath.calibdata}{dataNight}/S_{setNum:02d}/"
    midpointImage = f"{imgsetp}ib022.fit"
    obsTime, obsSite = get_site_info(midpointImage)

    # Load in catalg of SAO stars with V <= 7.5 mag
    saoStarsFile = f"{filepath.spreadsheets}SAO_stars75.xlsx"
    saoStars = pd.read_excel(saoStarsFile)

    # Get SAO star coordinates
    saoCoord = coord.SkyCoord(
        ra=saoStars['RA_deg'].values,
        dec=saoStars['Dec_deg'].values,
        unit=(u.deg, u.deg),
        frame='icrs'
    )

    # Convert to Alt-Az topocentric coordinates
    saoTopoCoord = saoCoord.transform_to(
        coord.AltAz(obstime=obsTime, location=obsSite)
    )
    saoStars['altitude'] = saoTopoCoord.alt.deg

    # Shorten list to stars above horizon
    saoStarsOnSky = saoStars.loc[
        saoTopoCoord.alt.deg > 0
    ].reset_index(drop=True)

    # Get Earth-Sun distance in Astronomical Units
    sunCoord = coord.get_body('sun', obsTime, obsSite)
    distEarthToSun = sunCoord.distance.value

    # Get coordinate and brightness info for planets
    planets = [] # empty list to store planet info
    planetNames = ['venus','mars','jupiter','saturn']
    for i,p in enumerate(planetNames):

        # Get planet apparent magnitude and Alt-Az coordinates
        mag, az, alt = planet_brightness(
            p, obsTime, obsSite, distEarthToSun
        )

        # Generate dataframe entry
        pdf = pd.DataFrame({
            'planet': p,
            'azimuth': az,
            'altitude': alt,
            'Mag': mag
        },index=[i])
        planets.append(pdf)

    # Concatenate into single dataframe and remove planets below horizon
    planets = pd.concat(planets)
    planetsOnSky = planets.loc[
        planets.altitude > 0
    ].reset_index(drop=True)

    # Calculate airmass for stars and planets
    saoStarsOnSky['ZA'] = 90 - saoStarsOnSky.altitude
    planetsOnSky['ZA'] = 90 - planetsOnSky.altitude
    sAirmass = calculate_airmass(saoStarsOnSky.ZA.values)
    pAirmass = calculate_airmass(planetsOnSky.ZA.values)
    saoStarsOnSky['airmass'] = sAirmass
    planetsOnSky['airmass'] = pAirmass

    # Calculate extincted magnitudes
    saoStarsOnSky['extMag'] = saoStarsOnSky.Vmag + extCoeff * sAirmass
    planetsOnSky['extMag'] = planetsOnSky.Mag + extCoeff * pAirmass

    # Calculate net flux to SQM, setting values to 0 at Zenith Angles
    costermStars = n.cos(n.deg2rad(saoStarsOnSky.ZA))
    costermPlanets = n.cos(n.deg2rad(planetsOnSky.ZA))
    saoStarsOnSky['netflux_sqm'] = costermStars * 10**(-saoStarsOnSky.extMag/2.5)
    planetsOnSky['netflux_sqm'] = costermPlanets * 10**(-planetsOnSky.extMag/2.5)
    saoStarsOnSky.loc[saoStarsOnSky.ZA > 60, 'netflux_sqm'] = 0.0
    planetsOnSky.loc[planetsOnSky.ZA > 60, 'netflux_sqm'] = 0.0

    # Calculate total flux values
    totalFluxStars = saoStarsOnSky.netflux_sqm.sum()
    totalFluxPlanets = planetsOnSky.netflux_sqm.sum()

    # Calculate star and planet flux per square arcsecond
    deg2sec = 3600*3600 # square arcseconds per square degree
    arcsecOnSky = 2*n.pi*deg2sec*57.296*(57.296-57.3/n.sqrt(3))
    starFlux = totalFluxStars / arcsecOnSky
    planetFlux = totalFluxPlanets / arcsecOnSky

    # Calculate synthetic SQM
    backgroundFlux = 10**(-skySQM/2.5)
    sqm = -2.5 * n.log10(starFlux + planetFlux + backgroundFlux)

    return sqm


def append_to_calibreport(metrics):
    '''
    Functin to add sky brightness metrics to
    calibreport Excel file.

    Parameters:
    -----------
    metrics: DataFrame
        DataFrame containing calculated metrics
    '''

    # Get path to calibreport Excel file
    dataNight = metrics['datanight'].iloc[0]
    excelFile = f"{filepath.calibdata}{dataNight}/calibreport.xlsx"

    # Open calibreport for writing
    with pd.ExcelWriter(
        excelFile, 
        engine='openpyxl', 
        if_sheet_exists='overlay', 
        mode='a') as writer:

        # Grab the Report sheet
        worksheet = writer.sheets['Report']

        # Iterate over each data set
        for i,row in metrics.iterrows():

            # Skip over B-band rows
            if row['filter'] == 'B':
                continue

            # Get the set number
            setnum = row['dataset']

            # Add metrics to Excel sheet
            worksheet.cell(row=21+2*(setnum-1), column=11, value=row['zenith_mag'])
            worksheet.cell(row=21+2*(setnum-1), column=12, value=row['brightest_mag'])
            worksheet.cell(row=21+2*(setnum-1), column=13, value=row['faintest_mag'])
            worksheet.cell(row=21+2*(setnum-1), column=14, value=row['SQM_synthetic'])


#------------------------------------------------------------------------------#
#-------------------              Main Program              -------------------#
#------------------------------------------------------------------------------#

def calculate_sky_quality(dnight,sets,filter,albedo):
    '''
    Main program for computing SQI and SQM metrics and simple
    aperture photometry metrics.

    Parameters:
    -----------
    dnight: str
        Name of data night to process
    sets: list
        List of data sets to process
    filter: str
        Name of photometric filter
    albedo: float
        Site-specific albedo

    Returns:
    --------
    sqOutput: DataFrame
        Dataframe containing light pollution metrics
    '''

    # Filter paths
    F = {'V':'', 'B':'B/'}

    # Load in zeropoint, platescale, and exposure time values
    calsetp = f"{filepath.calibdata}{dnight}/"
    extinctionFile = f"{calsetp}extinction_fit_{filter}.xlsx"
    extinctionData  = pd.read_excel(extinctionFile)
    zeropoint = extinctionData['zeropoint_default'].values
    platescale = extinctionData['avg_scale'].mean()
    exptime = extinctionData['exptime'].values

    # Calculate image scale offset
    scaleOffset = 2.5*n.log10((platescale*60)**2)

    # Loop through each data set
    sqOutput = []
    for s in sets:

        # Set path for grid datasets
        setnum = int(s[0])
        gridsetp = f"{filepath.griddata}{dnight}/S_{setnum:02d}/{F[filter]}"

        # Print status update
        print(f"{PREFIX}Processing {dnight} Set-{setnum}")

        # Get zeropoint and exposure time for given dataset
        zpSet = zeropoint[setnum-1]
        expSet = exptime[setnum-1]

        # Calculate sky quality indices for each mask
        sqiAllsky = calc_SQI(gridsetp, 'horizon')
        sqiZ80 = calc_SQI(gridsetp, 'ZA80')
        sqiZ70 = calc_SQI(gridsetp, 'ZA70')

        # Calculate sky background photometry and synthetic SQMs
        skySQM, imgPhot = calc_skyonly_SQM(dnight, setnum, filter)
        synSQM = calc_synthetic_SQM(dnight, setnum, filter, skySQM)

        # Measure brightness at center of zenith images
        zenithPhot = measure_zenith_images(f"{calsetp}/S_{setnum:02d}/")
        zenithPhot['mag'] = zpSet + scaleOffset - zenithPhot['inst_mag']
        magZenith2 = zenithPhot[zenithPhot.image == 'zenith2'].mag.iloc[0]

        # Measure All-sky and above ZA-70 total magnitudes
        aduSum = imgPhot.ADU.sum()
        aduSumZA70 = imgPhot[imgPhot.Pany >= 20].ADU.sum()
        scaleFactor = 4 * 3600 * 3600 / (platescale*60)**2
        allskyMag = zpSet - 2.5*n.log10(aduSum * scaleFactor / expSet)
        za70Mag = zpSet - 2.5*n.log10(aduSumZA70 * scaleFactor / expSet)

        # Get brightest and faintest magnitudes measured
        imgPhot['mag'] = zpSet + scaleOffset - 2.5*n.log10(imgPhot['ADU'] / expSet)
        faintestMag = imgPhot[imgPhot.Pany >= 46].mag.quantile(0.98)
        brightestMag = imgPhot.mag.min()

        # Calculate scalar illuminance
        scalarIllum = 0.00254 * 10**(-0.4*allskyMag) * albedo / 4

        # Print out results
        print(f'{PREFIX}SQI to Observed Horizon = {sqiAllsky:.2f}')
        print(f'{PREFIX}SQI to Zenith Angle 80  = {sqiZ80:.2f}')
        print(f'{PREFIX}SQI to Zenith Angle 70  = {sqiZ70:.2f}')
        print(f'{PREFIX}Sky Background SQM      = {skySQM:.2f}')
        print(f'{PREFIX}Full Synthetic SQM      = {synSQM:.2f}')

        # Generate dataframe entry for given dataset
        sqEntry = pd.DataFrame(
            {
                'datanight': dnight,
                'dataset': setnum,
                'filter': filter,
                'SQI_allsky': sqiAllsky,
                'SQI_ZA80': sqiZ80,
                'SQI_ZA70': sqiZ70,
                'SQM_sky': skySQM,
                'SQM_synthetic': synSQM,
                'zenith_mag': magZenith2,
                'allsky_mag': allskyMag,
                'za70_mag': za70Mag,
                'brightest_mag': brightestMag,
                'faintest_mag': faintestMag,
                'scalar_illum': scalarIllum
            },
            index = [setnum-1]
        )
        sqOutput.append(sqEntry)

    # Create final dataframe output
    sqOutput = pd.concat(sqOutput)

    # Add metrics to calibreport
    append_to_calibreport(sqOutput)

    return sqOutput
