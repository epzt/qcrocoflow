import numpy as np
import datetime
import netCDF4 as nc
import os
from .ap2ep_Script import ap2ep
from .interpTPXO_Script import interpTPXO
from .readtpxodata_Script import readTPXOdata
from .TPXOnodalfactors_Script import TPXOnodalfactors
from .Vphase_Script import Vphase
from datetime import date , datetime
from PyQt5.QtWidgets import QMessageBox




def TPXO2CROCO(t0, ROMSnames, crocoGRD, fnOut, ndays, TPXO_path, grdname):
    """
    Convert TPXO data to forcing Croco format.
        t0: Initial date for tide time, in MATLAB datenum format.
        ROMSnames: List of harmonic components.
        crocoGRD: Path to the ROMS grid file.
        fnOut: Path to the output file.
        ndays: Number of days for which data is to be extracted.
        TPXO_path: Path to the TPXO data directory.
    """

    # Convert the MATLAB date to a Python datetime object
    date_1968 = datetime(1968, 5, 23)
    datenum_1968 = date_1968.toordinal() + 366
    python_datetime = datetime.fromordinal(int(t0) - 366)
    ROMStitle = f'ROMS TPXO data for {python_datetime}'

    # Open the croco grid file and extract relevant variables
    with nc.Dataset(crocoGRD, 'r') as grid:
        lonR = np.mod(grid.variables['lon_rho'][:], 360)
        latR = grid.variables['lat_rho'][:]
        maskR = grid.variables['mask_rho'][:]
        maskP = grid.variables['mask_psi'][:]

    L, M = np.shape(maskP)

    # Define the dictionary for the TPXO data paths
    # Each harmonic component has a separate file for elevation and velocity data
    TPXOfile = {
        'grid': {
            'atlas30': os.path.join(TPXO_path, 'grid_tpxo9_atlas_30_v5.nc'),
        },
        'elev': {
            'k1': os.path.join(TPXO_path, 'hf.k1_tpxo9_atlas_30_v5.nc'),
            'k2': os.path.join(TPXO_path, 'hf.k2_tpxo9_atlas_30_v5.nc'),
            'm2': os.path.join(TPXO_path, 'hf.m2_tpxo9_atlas_30_v5.nc'),
            'm4': os.path.join(TPXO_path, 'hf.m4_tpxo9_atlas_30_v5.nc'),
            'mf': os.path.join(TPXO_path, 'hf.mf_tpxo9_atlas_30_v5.nc'),
            'mm': os.path.join(TPXO_path, 'hf.mm_tpxo9_atlas_30_v5.nc'),
            'mn4': os.path.join(TPXO_path, 'hf.mn4_tpxo9_atlas_30_v5.nc'),
            'ms4': os.path.join(TPXO_path, 'hf.ms4_tpxo9_atlas_30_v5.nc'),
            'n2': os.path.join(TPXO_path, 'hf.n2_tpxo9_atlas_30_v5.nc'),
            'o1': os.path.join(TPXO_path, 'hf.o1_tpxo9_atlas_30_v5.nc'),
            'p1': os.path.join(TPXO_path, 'hf.p1_tpxo9_atlas_30_v5.nc'),
            'q1': os.path.join(TPXO_path, 'hf.q1_tpxo9_atlas_30_v5.nc'),
            's2': os.path.join(TPXO_path, 'hf.s2_tpxo9_atlas_30_v5.nc'),
        },
        'vel': {
            'k1': os.path.join(TPXO_path, 'uv.k1_tpxo9_atlas_30_v5.nc'),
            'k2': os.path.join(TPXO_path, 'uv.k2_tpxo9_atlas_30_v5.nc'),
            'm2': os.path.join(TPXO_path, 'uv.m2_tpxo9_atlas_30_v5.nc'),
            'm4': os.path.join(TPXO_path, 'uv.m4_tpxo9_atlas_30_v5.nc'),
            'mf': os.path.join(TPXO_path, 'uv.mf_tpxo9_atlas_30_v5.nc'),
            'mm': os.path.join(TPXO_path, 'uv.mm_tpxo9_atlas_30_v5.nc'),
            'mn4': os.path.join(TPXO_path, 'uv.mn4_tpxo9_atlas_30_v5.nc'),
            'ms4': os.path.join(TPXO_path, 'uv.ms4_tpxo9_atlas_30_v5.nc'),
            'n2': os.path.join(TPXO_path, 'uv.n2_tpxo9_atlas_30_v5.nc'),
            'o1': os.path.join(TPXO_path, 'uv.o1_tpxo9_atlas_30_v5.nc'),
            'p1': os.path.join(TPXO_path, 'uv.p1_tpxo9_atlas_30_v5.nc'),
            'q1': os.path.join(TPXO_path, 'uv.q1_tpxo9_atlas_30_v5.nc'),
            's2': os.path.join(TPXO_path, 'uv.s2_tpxo9_atlas_30_v5.nc'),
        }
    }

    print('Reading TPXO')
    TPXO = readTPXOdata(TPXOfile, ROMSnames, lonR, latR)

    # Define the dictionary for the TPXO harmonics metadata
    # Each harmonic has a set of Doodson numbers, reference phase, and speed
    TPXOharmonics = {
        'mm': [0, 1, 0, -1, 0, 0, 0, 0.5443747],
        'mf': [0, 2, 0, 0, 0, 0, 0, 1.0980331],
        'q1': [1, -3, 1, 1, 0, 0, 270, 13.3986607],
        'o1': [1, -2, 1, 0, 0, 0, 270, 13.9430351],
        'p1': [1, 0, -1, 0, 0, 0, 270, 14.9589310],
        'k1': [1, 0, 1, 0, 0, 0, 90, 15.0410690],
        'n2': [2, -3, 2, 1, 0, 0, 0, 28.4397297],
        'm2': [2, -2, 2, 0, 0, 0, 0, 28.9841042],
        's2': [2, 0, 0, 0, 0, 0, 0, 30.0000000],
        'k2': [2, 0, 2, 0, 0, 0, 0, 30.0821381],
        'mn4': [4, -5, 4, 1, 0, 0, 0, 57.4238319],
        'm4': [4, -4, 4, 0, 0, 0, 0, 57.9682083],
        'ms4': [4, -2, 2, 0, 0, 0, 0, 58.9841042],
    }



    # Initialize dictionaries and lists to store harmonic information
    ROMSharmonics = {}
    Vdeg = {}
    Nharmonics = len(ROMSnames)
    ROMSperiods = [0] * Nharmonics
    for n in range(Nharmonics):
        ROMSharmonics[ROMSnames[n]] = TPXOharmonics[ROMSnames[n]]
        ROMSperiods[n] = 360 / TPXOharmonics[ROMSnames[n]][-1]

        Vdeg[ROMSnames[n]] = Vphase(t0, ROMSharmonics[ROMSnames[n]])

    # Compute the nodal corrections for the chosen date
    fFac, uFac = TPXOnodalfactors(t0 + ndays / 2, ROMSnames)

    # Initialize arrays to store TPXO data
    zamp, zpha, uamp, upha, vamp, vpha, major, eccentricity, inclination, phase = (
        np.zeros((L + 1, M + 1, Nharmonics)).transpose(1, 0, 2) for _ in range(10)
    )

    # Loop over each harmonic and interpolate the TPXO data onto the ROMS grid
    count_a30 = -1
    for k, name in enumerate(ROMSnames):
        count_a30 += 1


        harmonic = TPXO['harmonic_a30'][count_a30]['harmonic']
        print(f'Interpolating {harmonic} amplitudes')

        ei = interpTPXO(TPXO['h'], count_a30, lonR, latR, maskR)
        zamp[:, :, k] = np.abs(ei) * fFac[name]
        zpha[:, :, k] = np.mod(-np.angle(ei) * 180 / np.pi - uFac[name] - Vdeg[name], 360)

        print(f'Interpolating {harmonic} U')
        ei = interpTPXO(TPXO['U'], count_a30, lonR, latR, maskR)
        uamp[:, :, k] = np.abs(ei) * fFac[name]
        upha[:, :, k] = np.mod(-np.angle(ei) * 180 / np.pi - uFac[name] - Vdeg[name], 360)

        print(f'Interpolating {harmonic} V')
        ei = interpTPXO(TPXO['V'], count_a30, lonR, latR, maskR)
        vamp[:, :, k] = np.abs(ei) * fFac[name]
        vpha[:, :, k] = np.mod(-np.angle(ei) * 180 / np.pi - uFac[name] - Vdeg[name], 360)

        maj, ecc, inc, pha = ap2ep(np.squeeze(uamp[:, :, k]), np.squeeze(upha[:, :, k]), np.squeeze(vamp[:, :, k]),
                                   np.squeeze(vpha[:, :, k]))
        ecc[np.isnan(ecc)] = 0
        major[:, :, k] = maj
        eccentricity[:, :, k] = ecc
        inclination[:, :, k] = inc
        phase[:, :, k] = pha

    # Compute the minor axis from the major axis and eccentricity
    minor = major * eccentricity


    # Set up output netCDF forcing file
    print('Creating netcdf file')

    # global attributes
    varname = ['tide_period', 'tide_Ephase', 'tide_Eamp', 'tide_Cmin', 'tide_Cmax',
               'tide_Cangle', 'tide_Cphase', 'tide_Uamp', 'tide_Uphase',
               'tide_Vamp', 'tide_Vphase']
    vartype = ['f8'] * len(varname)
    vardims = [('tide_period',)] + [('tide_period','eta_rho', 'xi_rho')] * 10
    nVar = len(varname)

    # Create netCDF4 Dataset
    ds = nc.Dataset(fnOut, 'w', format='NETCDF4')

    # Define dimensions
    ds.createDimension('tide_period', Nharmonics)
    ds.createDimension('xi_rho', M + 1)
    ds.createDimension('eta_rho', L + 1)

    # Define variables
    for i in range(nVar):
        ds.createVariable(varname[i], vartype[i], vardims[i])

    ROMSnames_att = " ".join(ROMSnames)
    today = date.today().strftime("%Y%m%d")

    ds.title = ROMStitle
    ds.Creation_date = today
    ds.grd_file = crocoGRD
    ds.type = 'ROMS forcing file from TPXO'
    ds.ini_date_datenumber = t0
    ds.ini_date_mjd = (t0 - datenum_1968)
    ds.components = ROMSnames_att

    # Assign attributes to each variable
    attnames = [['long_name', 'units']] * nVar
    attvals = [['Tide angular period', 'hours'],
               ['Tide elevation phase angle', 'degrees'],
               ['Tide elevation amplitude', 'meters'],
               ['Tidal current ellipse semi-minor axis', 'meter second-1'],
               ['Tidal current ellipse semi-major axis', 'meter second-1'],
               ['Tidal current ellipse inclination angle', 'degrees between semi-major axis and east'],
               ['Tidal current phase angle', 'degrees'],
               ['Tidal current U-component amplitude', 'meters'],
               ['Tidal current U-component phase', 'degrees'],
               ['Tidal current V-component amplitude', 'meters'],
               ['Tidal current V-component phase', 'degrees']]

    for n in range(nVar):
        for i in range(len(attnames[n])):
            ds[varname[n]].setncattr(attnames[n][i], attvals[n][i])

    # Write data to the output file
    # We transpose our data arrays to match the dimensions in our netCDF dataset
    ds['tide_period'][:] = ROMSperiods
    ds['tide_Eamp'][:] = np.transpose(zamp, (2, 1, 0))
    ds['tide_Ephase'][:] = np.transpose(zpha, (2, 1, 0))
    ds['tide_Cmax'][:] = np.transpose(major, (2, 1, 0))
    ds['tide_Cmin'][:] = np.transpose(minor, (2, 1, 0))
    ds['tide_Cangle'][:] = np.transpose(inclination, (2, 1, 0))
    ds['tide_Cphase'][:] = np.transpose(phase, (2, 1, 0))
    ds['tide_Uamp'][:] = np.transpose(uamp, (2, 1, 0))
    ds['tide_Uphase'][:] = np.transpose(upha, (2, 1, 0))
    ds['tide_Vamp'][:] = np.transpose(vamp, (2, 1, 0))
    ds['tide_Vphase'][:] = np.transpose(vpha, (2, 1, 0))

    # Close the netCDF dataset
    ds.close()

    QMessageBox.information(None, "Interpolation Completed", f"Interpolation completed, the forcing file was successfully generated with the name {grdname}_tides.nc")





