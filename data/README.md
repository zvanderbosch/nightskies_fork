# Data (NSNSD CCD Processing Pipeline)

## [CCD.zip](./CCD.zip)

This ZIP file contains the directory structure assumed by the pipeline. The directory structure is shown below along with a description of the primary contents of each folder to the right:
```
CCD 
└─── Data                          -Contains the filelist.txt processing list
│    └─── calibdata                -Calibrated images and data validation outputs
│    └─── fielddata                -Raw images and calibration files
│    └─── graphics                 -Final all-sky skyglow images per data set
│    └─── griddata                 -ArcGIS grids and mosaics per data set
│    └─── maps                     -ArcGIS map templates for producing graphics
│    └─── rasters                  -Master ArcGIS grids for natural sky modeling
│    │    └─── scratch_fullres     -ArcGIS workspace for full-resolution mosaics
│    │    └─── scratch_galactic    -ArcGIS workspace for Galactic model mosaics
│    │    └─── scratch_median      -ArcGIS workspace for median-filtered mosaics
│    │    └─── scratch_metrics     -ArcGIS workspace for claculating skyglow metrics
│    │    └─── scratch_natsky      -ArcGIS workspace for natural sky mosaics
│    │    └─── scratch_zodiacal    -ArcGIS workspace for Zodiacal model mosaics
│    └─── standards                -Standard star catalogs for photometric calibration
│    └─── tables                   -Summary tables with metadata and light pollution metrics
│ 
└─── Images
     └─── Linearity Curves         -Master linearity curves per CCD camera
     └─── Master                   -Master Flat/Bias/Thermal images per CCD camera
```

## [filelist.txt](./filelist.txt)

This file tells the pipeline which data sets to process. It should be placed in the `CCD --> Data` directory and has the following fields:
   - `Dataset`: Name of data night to process (e.g. ROMO241004)
   - `V_band`: Yes or No, whether to process V-band images
   - `B_band`: Yes or No, whether to process B-band images
   - `Flat_V`: Name of master flat file used to calibrate V-band images
   - `Flat_B`: Name of master flat file used to calibrate B-band images
   - `Curve`: Name of linearity response curve file used to calibrated images
   - `Zeropoint`: The default zeropoint (mag) for the CCD camera used
   - `Processor`: Name of data processor with an **underscore** between first initial and last name (e.g. J_Doe)
   - `Central_AZ`: Azimuth coordinate to place at the center of final panoramic graphics
   - `Location`: Descriptive park name (e.g. Rocky_Mountain_NP), using **underscores** instead of spaces
