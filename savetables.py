#-----------------------------------------------------------------------------#
#savetables.py
#
#NPS Night Skies Program
#
#Last updated: 2025/05/28
#
#This script generates the final output tables consolidating
#data night attributes, calibration, and light pollution 
#metrics for each data set.
#
#Note: 
#
#Input:
#   (1) 
#
#Output:
#   (1) 
#
#History:
#	Zach Vanderbosch -- Created script
#
#-----------------------------------------------------------------------------#

import os
import stat
import numpy as n
import pandas as pd
import matplotlib.pyplot as plt

# Local Source
import filepath
import printcolors as pc

#------------------------------------------------------------------------------#
#-------------------            Define Functions            -------------------#
#------------------------------------------------------------------------------#




#------------------------------------------------------------------------------#
#-------------------              Main Program              -------------------#
#------------------------------------------------------------------------------#

def generate_tables(dnight):

    excelFile = f"{filepath.tables}{dnight}.xlsx"

    # Initial creation of excel file
    with pd.ExcelWriter(excelFile, engine='openpyxl', mode='w') as writer:

        # Create an openpyxl workbook object
        workbook = writer.book

        # Add each sheet
        workbook.create_sheet("NIGHT METADATA")
        workbook.create_sheet("CITIES")
        workbook.create_sheet("SET METADATA")
        workbook.create_sheet("CALIBRATION")
        workbook.create_sheet("EXTINCTION")
        workbook.create_sheet("IMG COORDS")
        workbook.create_sheet("V2 PHOTOMETRY")
        workbook.create_sheet("V4 PHOTOMETRY")
        workbook.create_sheet("GLARE")
        workbook.create_sheet("NATSKY")
        workbook.create_sheet("LP")
        workbook.create_sheet("LP80")
        workbook.create_sheet("LP70")
        workbook.create_sheet("ZONES")
        workbook.create_sheet("V4 PERCENTILES ALL")
        workbook.create_sheet("V4 PERCENTILES LP")