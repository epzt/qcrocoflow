import numpy as np
import xarray as xr
import os
from .Qsat import qsat
from netCDF4 import Dataset


#------------------------------------------------------------------#
def interp_ERA5(ATMO_dir, Y, M, nc_blk, angle, tlen):

    # Air temperature
    vname='T2M'
    with xr.open_dataset(os.path.join(ATMO_dir, f"{vname}_Y{Y}M{M}.nc")) as ds:
        tair = ds[vname].mean(dim=['lon', 'lat']) - 273.15

    # Relative humidity
    vname='Q'
    with xr.open_dataset(os.path.join(ATMO_dir, f"{vname}_Y{Y}M{M}.nc")) as ds:
        shum = ds[vname].mean(dim=['lon', 'lat'])
    rhum = shum / qsat(tair)  # Assuming qsat() function is defined

    # Precipitation rate
    vname='TP'
    with xr.open_dataset(os.path.join(ATMO_dir, f"{vname}_Y{Y}M{M}.nc")) as ds:
        prate = ds[vname].mean(dim=['lon', 'lat']) * 0.1 * (24.0 * 60.0 * 60.0)
    prate = np.where(prate < 1.e-4, 0, prate)

    # Shortwave flux
    vname='SSR'
    with xr.open_dataset(os.path.join(ATMO_dir, f"{vname}_Y{Y}M{M}.nc")) as ds:
        radsw = ds[vname].mean(dim=['lon', 'lat'])
    radsw = np.where(radsw < 1.e-10, 0, radsw)

    # Longwave flux
    vname='STRD'
    with xr.open_dataset(os.path.join(ATMO_dir, f"{vname}_Y{Y}M{M}.nc")) as ds:
        dlwrf_in = ds[vname].mean(dim=['lon', 'lat'])

    # Wind
    vname='U10M'
    with xr.open_dataset(os.path.join(ATMO_dir, f"{vname}_Y{Y}M{M}.nc")) as ds:
        uwnd = ds[vname].mean(dim=['lon', 'lat'])

    vname='V10M'
    with xr.open_dataset(os.path.join(ATMO_dir, f"{vname}_Y{Y}M{M}.nc")) as ds:
        vwnd = ds[vname].mean(dim=['lon', 'lat'])

    wspd = np.sqrt(uwnd**2 + vwnd**2)

    # Rotations on the CROCO grid
    cosa = np.cos(angle)
    sina = np.sin(angle)
    u10 = (uwnd * cosa + vwnd * sina)
    u10 = u10.values.flatten()


    v10 = (vwnd * cosa - uwnd * sina)
    v10 = v10.values.flatten()


    # Fill the CROCO files

    nc_blk = Dataset(nc_blk,'a')
    for i in range(tlen):
        nc_blk['tair'][i, :, :] = tair[i]
        nc_blk['rhum'][i, :, :] = rhum[i]
        nc_blk['prate'][i, :, :] = prate[i]
        nc_blk['wspd'][i, :, :] = wspd[i]
        nc_blk['radlw_in'][i, :, :] = dlwrf_in[i]
        nc_blk['radsw'][i, :, :] = radsw[i]
        nc_blk['uwnd'][i, :, :] = u10[i]
        nc_blk['vwnd'][i, :, :] = v10[i]

