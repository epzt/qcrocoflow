import numpy as np
import netCDF4 as nc_module
import matplotlib.pyplot as plt

from netCDF4 import Dataset
from .create_grid import create_grid_file
from .get_metrics import get_metrics_file
from .add_topo import add_topo_file
from .process_mask import process_mask_file
from .uvp_mask import uvp_mask_file
from .get_angle import get_angle_file
from .smooth_grid import smoothgrid_file
import os

def make_grid_function(title, grdname, lon_min, lon_max, lat_min, lat_max, dl, hmin, topofile, directory):

    def rho2uvp(rfield):
        Mp, Lp = rfield.shape
        M = Mp - 1
        L = Lp - 1

        vfield = 0.5 * (rfield[0:M, :] + rfield[1:Mp, :])
        ufield = 0.5 * (rfield[:, 0:L] + rfield[:, 1:Lp])
        pfield = 0.5 * (ufield[0:M, :] + ufield[1:Mp, :])

        return ufield, vfield, pfield



    rtarget = 0.15


    theta = np.arctan((lat_max - lat_min) / (lon_max - lon_min))

    x = np.arange(lon_min, lon_max, dl)
    y = np.arange(lat_min, lat_max, dl)

    #coordonnées de base de la grille avant rotation
    Lonr_base, Latr_base = np.meshgrid(x, y)


    if theta != 0:
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        lons_diff = Lonr_base - lon_min
        lats_diff = Latr_base - lat_min

        Lonr = lon_min + lons_diff * cos_theta - lats_diff * sin_theta
        Latr = lat_min + lons_diff * sin_theta + lats_diff * cos_theta
    else:
        Lonr = Lonr_base
        Latr = Latr_base

    Lonu, Lonv, Lonp = rho2uvp(Lonr)
    Latu, Latv, Latp = rho2uvp(Latr)


    M, L = Latp.shape


    create_grid_file(L, M, grdname, title, directory)


    with nc_module.Dataset(grdname, 'a') as nc:
        nc['lat_u'][:] = Latu
        nc['lon_u'][:] = Lonu
        nc['lat_v'][:] = Latv
        nc['lon_v'][:] = Lonv
        nc['lat_rho'][:] = Latr
        nc['lon_rho'][:] = Lonr
        nc['lat_psi'][:] = Latp
        nc['lon_psi'][:] = Lonp

    print(' ')
    print(' Compute the metrics...')
    pm, pn, dndx, dmde = get_metrics_file(grdname)  

    L = pm.shape[1]
    M = pm.shape[0]

    xr = np.zeros_like(pm)
    yr = np.zeros_like(pm)

    for i in range(L-1):
        xr[:, i + 1] = xr[:, i] + 2 / (pm[:, i + 1] + pm[:, i])
    for j in range(M-1):
        yr[j + 1, :] = yr[j, :] + 2 / (pn[j + 1, :] + pn[j, :])
    xu, xv, xp = rho2uvp(xr)  
    yu, yv, yp = rho2uvp(yr)  

    dx = 1 / pm
    dy = 1 / pn

    dxmax = np.max(dx / 1000)
    dxmin = np.min(dx / 1000)
    dymax = np.max(dy / 1000)
    dymin = np.min(dy / 1000)

    print(' ')
    print(f' Min dx={dxmin:.3f} km - Max dx={dxmax:.3f} km')
    print(f' Min dy={dymin:.3f} km - Max dy={dymax:.3f} km')

    #get angle

    angle=get_angle_file(Latu,Lonu);
    # Coriolis parameter
    f = 4 * np.pi * np.sin(np.pi * Latr / 180) * 366.25 / (24 * 3600 * 365.25)

    # Fill the grid file
    print('Fill the grid file...')
    with Dataset(grdname, 'a') as nc:
        nc.variables['pm'][:] = pm
        nc.variables['pn'][:] = pn
        nc.variables['dndx'][:] = dndx
        nc.variables['dmde'][:] = dmde
        nc.variables['x_u'][:] = xu
        nc.variables['y_u'][:] = yu
        nc.variables['x_v'][:] = xv
        nc.variables['y_v'][:] = yv
        nc.variables['x_rho'][:] = xr
        nc.variables['y_rho'][:] = yr
        nc.variables['x_psi'][:] = xp
        nc.variables['y_psi'][:] = yp
        nc.variables['f'][:] = f
        nc.variables['angle'][:] = angle

    #add topo
    print("\nAdd topography...")
    grdname_full_path = os.path.join(directory, grdname)
    grdname_nc = grdname
    topofile_nc = topofile

    h = add_topo_file(grdname_nc, topofile_nc)
    plt.figure()
    plt.imshow(h, origin='lower')
    plt.title('Bathymetry')
    plt.colorbar()
    plt.show()



