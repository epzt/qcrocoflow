import netCDF4 as nc
import numpy as np
from datetime import datetime
import os
def create_grid_file(L, M, grdname, title, directory):
    Lp = L + 1
    Mp = M + 1
    
    # Création du chemin complet pour le fichier
    grdname_full_path = os.path.join(directory, grdname)
    nw = nc.Dataset(grdname, 'w', format='NETCDF4')

    # Create dimensions
    nw.createDimension('xi_u', L)
    nw.createDimension('eta_u', Mp)
    nw.createDimension('xi_v', Lp)
    nw.createDimension('eta_v', M)
    nw.createDimension('xi_rho', Lp)
    nw.createDimension('eta_rho', Mp)
    nw.createDimension('xi_psi', L)
    nw.createDimension('eta_psi', M)
    nw.createDimension('one', 1)
    nw.createDimension('two', 2)
    nw.createDimension('four', 4)
    nw.createDimension('bath', 1)

    # Create variables and attributes
    def create_ncvar(name, dimensions, long_name, units, dtype=np.double):
        var = nw.createVariable(name, dtype, dimensions)
        var.long_name = long_name
        var.units = units
        return var

    xl = create_ncvar('xl', 'one', 'domain length in the XI-direction', 'meter')
    el = create_ncvar('el', 'one', 'domain length in the ETA-direction', 'meter')
    depthmin = create_ncvar('depthmin', 'one', 'Shallow bathymetry clipping depth', 'meter')
    depthmax = create_ncvar('depthmax', 'one', 'Deep bathymetry clipping depth', 'meter')

    spherical = nw.createVariable('spherical', str, 'one')
    spherical.long_name = 'Grid type logical switch'
    spherical.option_T = 'spherical'

    angle = create_ncvar('angle', ('eta_rho', 'xi_rho'), 'angle between xi axis and east', 'radian')
    h = create_ncvar('h', ('eta_rho', 'xi_rho'), 'Final bathymetry at RHO-points', 'meter')
    hraw = create_ncvar('hraw', ('bath', 'eta_rho', 'xi_rho'), 'Working bathymetry at RHO-points', 'meter')
    alpha = create_ncvar('alpha', ('eta_rho', 'xi_rho'), 'Weights between coarse and fine grids at RHO-points', '')
    f = create_ncvar('f', ('eta_rho', 'xi_rho'), 'Coriolis parameter at RHO-points', 'second-1')
    pm = create_ncvar('pm', ('eta_rho', 'xi_rho'), 'curvilinear coordinate metric in XI', 'meter-1')
    pn = create_ncvar('pn', ('eta_rho', 'xi_rho'), 'curvilinear coordinate metric in ETA', 'meter-1')
    dndx = create_ncvar('dndx', ('eta_rho', 'xi_rho'), 'xi derivative of inverse metric factor pn', 'meter')
    dmde = create_ncvar('dmde', ('eta_rho', 'xi_rho'), 'eta derivative of inverse metric factor pm', 'meter')

    nw.createVariable('x_rho', 'f8', ('eta_rho', 'xi_rho'))
    nw['x_rho'].long_name = 'x location of RHO-points'
    nw['x_rho'].units = 'meter'

    nw.createVariable('x_u', 'f8', ('eta_u', 'xi_u'))
    nw['x_u'].long_name = 'x location of U-points'
    nw['x_u'].units = 'meter'

    nw.createVariable('x_v', 'f8', ('eta_v', 'xi_v'))
    nw['x_v'].long_name = 'x location of V-points'
    nw['x_v'].units = 'meter'

    nw.createVariable('x_psi', 'f8', ('eta_psi', 'xi_psi'))
    nw['x_psi'].long_name = 'x location of PSI-points'
    nw['x_psi'].units = 'meter'

    # y location variables
    nw.createVariable('y_rho', 'f8', ('eta_rho', 'xi_rho'))
    nw['y_rho'].long_name = 'y location of RHO-points'
    nw['y_rho'].units = 'meter'

    nw.createVariable('y_u', 'f8', ('eta_u', 'xi_u'))
    nw['y_u'].long_name = 'y location of U-points'
    nw['y_u'].units = 'meter'

    nw.createVariable('y_v', 'f8', ('eta_v', 'xi_v'))
    nw['y_v'].long_name = 'y location of V-points'
    nw['y_v'].units = 'meter'

    nw.createVariable('y_psi', 'f8', ('eta_psi', 'xi_psi'))
    nw['y_psi'].long_name = 'y location of PSI-points'
    nw['y_psi'].units = 'meter'

    # longitude variables
    nw.createVariable('lon_rho', 'f8', ('eta_rho', 'xi_rho'))
    nw['lon_rho'].long_name = 'longitude of RHO-points'
    nw['lon_rho'].units = 'degree_east'

    nw.createVariable('lon_u', 'f8', ('eta_u', 'xi_u'))
    nw['lon_u'].long_name = 'longitude of U-points'
    nw['lon_u'].units = 'degree_east'

    nw.createVariable('lon_v', 'f8', ('eta_v', 'xi_v'))
    nw['lon_v'].long_name = 'longitude of V-points'
    nw['lon_v'].units = 'degree_east'

    nw.createVariable('lon_psi', 'f8', ('eta_psi', 'xi_psi'))
    nw['lon_psi'].long_name = 'longitude of PSI-points'
    nw['lon_psi'].units = 'degree_east'

    # latitude variables
    nw.createVariable('lat_rho', 'f8', ('eta_rho', 'xi_rho'))
    nw['lat_rho'].long_name = 'latitude of RHO-points'
    nw['lat_rho'].units = 'degree_north'

    nw.createVariable('lat_u', 'f8', ('eta_u', 'xi_u'))
    nw['lat_u'].long_name = 'latitude of U-points'
    nw['lat_u'].units = 'degree_north'

    nw.createVariable('lat_v', 'f8', ('eta_v', 'xi_v'))
    nw['lat_v'].long_name = 'latitude of V-points'
    nw['lat_v'].units = 'degree_north'

    nw.createVariable('lat_psi', 'f8', ('eta_psi', 'xi_psi'))
    nw['lat_psi'].long_name = 'latitude of PSI-points'
    nw['lat_psi'].units = 'degree_north'

    nw.createVariable('mask_rho', 'f8', ('eta_rho', 'xi_rho'))
    nw['mask_rho'].long_name = 'mask on RHO-points'
    nw['mask_rho'].option_0 = 'land'
    nw['mask_rho'].option_1 = 'water'

    nw.createVariable('mask_u', 'f8', ('eta_u', 'xi_u'))
    nw['mask_u'].long_name = 'mask on U-points'
    nw['mask_u'].option_0 = 'land'
    nw['mask_u'].option_1 = 'water'

    nw.createVariable('mask_v', 'f8', ('eta_v', 'xi_v'))
    nw['mask_v'].long_name = 'mask on V-points'
    nw['mask_v'].option_0 = 'land'
    nw['mask_v'].option_1 = 'water'

    nw.createVariable('mask_psi', 'f8', ('eta_psi', 'xi_psi'))
    nw['mask_psi'].long_name = 'mask on PSI-points'
    nw['mask_psi'].option_0 = 'land'
    nw['mask_psi'].option_1 = 'water'

    # Create global attributes
    nw.title = title
    date = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    nw.date = date
    nw.type = 'CROCO grid file'

    nw.close()