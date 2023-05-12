import numpy as np
import netCDF4 as nc
from .spheric_dist import spheric_dist_file

def get_metrics_file(grdname):
    # Read in the grid
    with nc.Dataset(grdname, 'r') as ds:
        latu = ds['lat_u'][:]
        lonu = ds['lon_u'][:]
        latv = ds['lat_v'][:]
        lonv = ds['lon_v'][:]

    Mp, L = latu.shape
    M, Lp = latv.shape
    Lm = L - 1
    Mm = M - 1

    # pm and pn
    dx = np.zeros((Mp, Lp))
    dy = np.zeros((Mp, Lp))
    dx[:, 1:L] = spheric_dist_file(latu[:, 0:Lm], latu[:, 1:L], lonu[:, 0:Lm], lonu[:, 1:L])
    dx[:, 0] = dx[:, 1]
    dx[:, Lp - 1] = dx[:, L - 1]

    dy[1:M, :] = spheric_dist_file(latv[0:Mm, :], latv[1:M, :], lonv[0:Mm, :], lonv[1:M, :])
    dy[0, :] = dy[1, :]
    dy[Mp - 1, :] = dy[M - 1, :]

    pm = 1.0 / dx
    pn = 1.0 / dy

    # dndx and dmde
    dndx = np.zeros((Mp, Lp))
    dmde = np.zeros((Mp, Lp))

    dndx[1:M, 1:L] = 0.5 * (1.0 / pn[1:M, 2:Lp] - 1.0 / pn[1:M, 0:Lm])
    dmde[1:M, 1:L] = 0.5 * (1.0 / pm[2:Mp, 1:L] - 1.0 / pm[0:Mm, 1:L])

    dndx[0, :] = 0
    dndx[Mp - 1, :] = 0
    dndx[:, 0] = 0
    dndx[:, Lp - 1] = 0

    dmde[0, :] = 0
    dmde[Mp - 1, :] = 0
    dmde[:, 0] = 0
    dmde[:, Lp - 1] = 0

    return pm, pn, dndx, dmde
