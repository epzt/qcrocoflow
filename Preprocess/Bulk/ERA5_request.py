#!/usr/bin/env python

# Script to download ECMWF ERA5 reanalysis datasets from the Climate Data
#  Store (CDS) of Copernicus https://cds.climate.copernicus.eu
#
#  This script use the CDS Phyton API[*] to connect and download specific ERA5
#  variables, for a chosen area and monthly date interval, required by CROCO to
#  perform simulations with atmospheric forcing. Furthermore, this script use
#  ERA5 parameter names and not parameter IDs as these did not result in stable
#  downloads.
#
#  Tested using Python 3.8.6 and Python 3.9.1. This script need the following
#  python libraries pre-installed: "calendar", "datetime", "json" and "os".
#
#  [*] https://cds.climate.copernicus.eu/api-how-to
#
#  Copyright (c) DDONOSO February 2021
#  e-mail:ddonoso@dgeo.udec.cl
#

#  You may see all available ERA5 variables at the following website
#  https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#ERA5:datadocumentation-Parameterlistings

# -------------------------------------------------
# Getting libraries and utilities
# -------------------------------------------------
import cdsapi
import os
import datetime
import calendar
from datetime import date
import json
from netCDF4 import Dataset as netcdf
import numpy as np
# -------------------------------------------------

from qgis.core import (QgsCoordinateReferenceSystem, QgsFeature, QgsGeometry,
                       QgsMapLayer, QgsPointXY, QgsProject, QgsRectangle,
                       QgsVectorLayer, QgsRasterLayer, QgsCoordinateTransform)
from qgis.utils import iface
from PyQt5.QtWidgets import QComboBox, QDialog, QFileDialog, QMessageBox, QLabel, QVBoxLayout, QApplication

# -------------------------------------------------


def ERA5_request_script(title, lon_min, lon_max, lat_min, lat_max, directory, Ymin, Ymax, Mmin, Mmax):
    transform_4326 = QgsCoordinateTransform(QgsCoordinateReferenceSystem('EPSG:3857'),
                                            QgsCoordinateReferenceSystem('EPSG:4326'),
                                            QgsProject.instance())
    QMessageBox.warning(None, "Warning",
                        "Warning: The download of ERA5 files may take a while. A message will notify you when the download is complete.")
    iface.actionShowPythonDialog().trigger()


    # Transforme les coordonnées lon_min, lon_max, lat_min, lat_max
    lon_min, lat_min = transform_4326.transform(lon_min, lat_min)
    lon_max, lat_max = transform_4326.transform(lon_max, lat_max)


    # -------------------------------------------------
    # Increment date by custom months
    # -------------------------------------------------
    def addmonths4date(date, addmonths):
        month = date.month - 1 + addmonths
        year = date.year + month // 12
        month = month % 12 + 1
        day = min(date.day, calendar.monthrange(year, month)[1])
        return datetime.date(year, month, day)

    ##config_dir = '../croco/Run_TEST/'           # must be the same than crocotools_param
    #
    config_name = title
    #
    # Original ERA5 directory
    #
    config_dir = os.path.join(directory, "BULK_FILE")
    era5_dir_raw = os.path.join(config_dir, "ERA5_native_" + config_name)
    era5_dir_processed = os.path.join(config_dir, "ERA5_" + config_name)
    #
    # extraction wave variables
    #
    wave_extract = False  # True to extract wave variables
    #
    # Dates limits
    #
    year_start = Ymin
    month_start = Mmin
    year_end = Ymax
    month_end = Mmax
    #
    # Year origin of time
    #
    Yorig = 2000
    #
    #
    # Request time (daily hours '00/01/.../23')
    #
    time = '00/01/02/03/04/05/06/07/08/09/10/11/12/13/14/15/16/17/18/19/20/21/22/23'
    #
    # Request variables (see available at ERA5_variables.json)
    variables = ['lsm', 'tp', 'strd', 'ssr', 't2m', 'q', 'u10', 'v10']  # note lsm is land_sea_mask
    #
    # Request area ([north, west, south, east])
    #
    ownArea = 1  # 0 if area from a crocotools_param.m file
    # 1 if own area

    lonmin = lon_min
    lonmax = lon_max
    latmin = lat_min
    latmax = lat_max
    #
    # Variable names and conversion coefficients
    # TP: convert from accumlated m in a hour into   kg m-2 s-1
    #
    cff_tp = 1000. / 3600.  # m in 1 hour -> kg m-2 s-1
    # Heat flux J m-2 in one hour into W m-2
    #
    cff_heat = 1. / 3600.  # J m-2 in 1 hour -> W m-2
    # Names, conversion coefficients and new units
    #
    variables = ['lsm', 'sst', 'tp', 'strd', 'ssr', 't2m', 'q', 'u10', 'v10']
    conv_cff = [1., 1., cff_tp, cff_heat, cff_heat, 1., 1., 1., 1.]
    units = ['(0-1)', 'K', 'kg m-2 s-1', 'W m-2', 'W m-2', 'K', 'kg kg-1', 'm s-1', 'm s-1']

    if wave_extract:
        ## append waves variables
        wave_var = ['swh', 'mwd', 'pp1d', 'cdww'];
        variables.extend(wave_var)
        wave_conv_cff = [1., 1., 1., 1.];
        conv_cff.extend(wave_conv_cff)
        wave_units = ['m', 'Degrees true', 's', 'dimensionless'];
        units.extend(wave_units)

    # *******************************************************************************
    #                         E N D     U S E R  *  O P T I O N S
    # *******************************************************************************

    # -------------------------------------------------
    dl = 0.5
    n_overlap = 0

    lonmin = str(float(lonmin) - dl)
    lonmax = str(float(lonmax) + dl)
    latmin = str(float(latmin) - dl)
    latmax = str(float(latmax) + dl)

    # -------------------------------------------------

    area = [latmax, lonmin, latmin, lonmax]

    # -------------------------------------------------
    # Setting raw output directory
    # -------------------------------------------------
    # Get the current directory
    os.makedirs(era5_dir_raw, exist_ok=True)

    # -------------------------------------------------
    # Loading ERA5 variables's information as
    # python Dictionary from JSON file
    # -------------------------------------------------
    script_dir = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(script_dir, 'ERA5_variables.json')

    with open(json_path, 'r') as jf:
        era5 = json.load(jf)

    # -------------------------------------------------
    # Downloading ERA5 datasets
    # -------------------------------------------------
    # Monthly dates limits
    monthly_date_start = datetime.datetime(year_start, month_start, 1)
    monthly_date_end = datetime.datetime(year_end, month_end, 1)

    # Length of monthly dates loop
    len_monthly_dates = (monthly_date_end.year - monthly_date_start.year) * 12 + \
                        (monthly_date_end.month - monthly_date_start.month) + 1

    # Initial monthly date
    monthly_date = monthly_date_start

    # Monthly dates loop
    for j in range(len_monthly_dates):

        # Year and month
        year = monthly_date.year;
        month = monthly_date.month;

        # Number of days in month
        days_in_month = calendar.monthrange(year, month)[1]

        # Date limits
        date_start = datetime.datetime(year, month, 1)
        date_end = datetime.datetime(year, month, days_in_month)
        # Ordinal date limits (days)
        n_start = datetime.date.toordinal(date_start)
        n_end = datetime.date.toordinal(date_end)
        # Overlapping date string limits (yyyy-mm-dd)
        datestr_start_overlap = datetime.date.fromordinal(n_start - n_overlap).strftime('%Y-%m-%d')
        datestr_end_overlap = datetime.date.fromordinal(n_end + n_overlap).strftime('%Y-%m-%d')
        # Overlapping date string interval
        vdate = datestr_start_overlap + '/' + datestr_end_overlap


        # Variables/Parameters loop
        for k in range(len(variables)):
            # Variable's name, long-name and level-type
            vname = variables[k]
            vlong = era5[vname][0]
            vlevt = era5[vname][3]

            # Request options
            options = {
                'product_type': 'reanalysis',
                'type': 'an',
                'date': vdate,
                'variable': vlong,
                'levtype': vlevt,
                'area': area,
                'format': 'netcdf',
            }

            # Add options to Variable without "diurnal variations"
            if vlong == 'sea_surface_temperature':
                options['time'] = '00'

            elif vlong == 'land_sea_mask':
                options['time'] = '00:00'

            else:
                options['time'] = time

            # Add options to Product "pressure-levels"
            if vlong == 'specific_humidity' or vlong == 'relative_humidity':
                options['pressure_level'] = '1000'
                product = 'reanalysis-era5-pressure-levels'

            # Product "single-levels"
            else:
                product = 'reanalysis-era5-single-levels'

            # Output filename
            fname = 'ERA5_ecmwf_' + vname.upper() + '_Y' + str(year) + 'M' + str(month).zfill(2) + '.nc'
            output = era5_dir_raw + '/' + fname

            # Information strings
            info_time_clock = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            info_monthly_date = monthly_date.strftime('%Y-%b')
            info_n_overlap = ' with ' + str(n_overlap) + ' overlapping day(s) '

            # Printing message on screen
            print('                                                           ')
            print('-----------------------------------------------------------')
            print('', info_time_clock, '                                        ')
            print(' Performing ERA5 data request, please wait...              ')
            print(' Date [yyyy-mmm] =', info_monthly_date + info_n_overlap, '   ')
            print(' Variable =', vlong, '                                       ')
            print('-----------------------------------------------------------')
            print('Request options: ')
            print(options)
            # Server ECMWF-API
            c = cdsapi.Client()
            # Do the request
            c.retrieve(product, options, output)
        # ---------------------------------------------------------------------
        # Next iteration to monthly date: add one month to current monthly date
        # ---------------------------------------------------------------------
        monthly_date = addmonths4date(monthly_date, 1)
    # !/usr/bin/env python
    #
    # ERA5_convert.py
    #
    # Script to convert ERA5 original files obtained from the Climate Data
    # Store (CDS) of Copernicus https://cds.climate.copernicus.eu into
    # a format and using units which can be used by the online interpolation of CROCO
    #
    # -------------------------------------------------
    # Setting processed output directory
    # -------------------------------------------------
    # Get the current directory
    os.makedirs(era5_dir_processed, exist_ok=True)

    # -------------------------------------------------
    # Loading ERA5 variables's information as
    # python Dictionary from JSON file
    # -------------------------------------------------

    with open(json_path, 'r') as jf:
        era5 = json.load(jf)
    #
    # -------------------------------------------------
    # Loop on Years and Months
    # -------------------------------------------------
    #
    for iyear in range(year_start, year_end + 1):
        for imonth in range(month_start, month_end + 1):
            #
            # -------------------------------------------------
            # Loop on variables names
            # -------------------------------------------------
            #
            for k in range(len(variables)):
                #
                # Variable's name, long-name and level-type
                #
                vname = variables[k]
                vlong = era5[vname][0]

                print('  Processing variable: ' + vname)
                #
                # Read input filedate.toordinal(date(Yorig,1,1))
                #

                fname_in = era5_dir_raw + '/ERA5_ecmwf_' + vname.upper() + '_Y' + str(iyear) + 'M' + str(imonth).zfill(
                    2) + '.nc'
                nc = netcdf(fname_in, 'r+', format='NETCDF4')
                time = nc.variables['time'][:]
                lat = nc.variables['latitude'][:]
                lon = nc.variables['longitude'][:]
                data = nc.variables[vname][:, :, :]
                nc.close()
                #
                # Flip latitudes (to have increasing latitudes...)
                #
                lat = np.flip(lat, axis=0)
                data = np.flip(data, axis=1)
                #
                # Missing values and multiply by cff to change unit
                #
                try:
                    mvalue = data.fill_value
                except AttributeError:
                    print('No fill value.. use nan')
                    mvalue = np.nan
                data = np.array(data)
                data = conv_cff[k] * data
                data[np.where(data == mvalue)] = 9999.

                #
                # Convert time from hours since 1900-1-1 0:0:0 into days since Yorig-1-1 0:0:0
                #

                time = time / 24.
                time = time - date.toordinal(date(Yorig, 1, 1)) \
                       + date.toordinal(date(1900, 1, 1))

                #
                # Changes names
                #

                if vname == 'u10':
                    vname = 'u10m'

                if vname == 'v10':
                    vname = 'v10m'
                #
                # Create and write output netcdf file
                #

                fname_out = era5_dir_processed + '/' + vname.upper() + '_Y' + str(iyear) + 'M' + str(imonth) + '.nc'

                nw = netcdf(fname_out, mode='w', format='NETCDF4')

                dimlon = nw.createDimension('lon', len(lon))
                dimlat = nw.createDimension('lat', len(lat))
                dimtime = nw.createDimension('time', None)

                varlon = nw.createVariable('lon', 'f4', ('lon',))
                varlat = nw.createVariable('lat', 'f4', ('lat',))
                vartime = nw.createVariable('time', 'f4', ('time',))
                vardata = nw.createVariable(vname.upper(), 'f4', ('time', 'lat', 'lon'))
                varlon.long_name = 'longitude of RHO-points'
                varlat.long_name = 'latitude of RHO-points'
                vartime.long_name = 'Time'
                varlon.units = 'degree_east'
                varlat.units = 'degree_north'
                vartime.units = 'days since ' + str(Yorig) + '-1-1'
                vardata.missing_value = 9999.
                vardata.units = units[k]
                vardata.long_name = vlong

                varlon[:] = lon
                varlat[:] = lat
                vartime[:] = time
                vardata[:] = data

                nw.close()
    #close python console
    iface.actionShowPythonDialog().trigger()

    # Print last message on screen
    QMessageBox.information(None, "Download Completed",
                            "Download completed successfully. You can now proceed with the creation of the bulk file.")


