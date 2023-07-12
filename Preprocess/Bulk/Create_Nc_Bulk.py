import netCDF4 as nc
import datetime
import os
def create_bulk(grdname, title, bulkt, bulkc, directory):
    grd = nc.Dataset(grdname, 'r')
    L = len(grd.dimensions['xi_psi'])
    M = len(grd.dimensions['eta_psi'])
    grd.close()
    Lp = L + 1
    Mp = M + 1
    frcname = directory

    nw = nc.Dataset(frcname, 'w', format='NETCDF4')
    nw.createDimension('xi_rho', Lp)
    nw.createDimension('eta_rho', Mp)
    nw.createDimension('xi_psi', L)
    nw.createDimension('eta_psi', M)
    nw.createDimension('xi_u', L)
    nw.createDimension('eta_u', Mp)
    nw.createDimension('xi_v', Lp)
    nw.createDimension('eta_v', M)
    nw.createDimension('bulk_time', None)

    var = nw.createVariable('bulk_time', 'f8', ('bulk_time',))
    var.long_name = 'bulk formulation execution time'
    var.units = 'days'
    var.cycle_length = bulkc

    # Create other variables similarly...
    var = nw.createVariable('tair', 'f8', ('bulk_time', 'eta_rho', 'xi_rho'))
    var.long_name = 'surface air temperature'
    var.units = 'Celsius'

    var = nw.createVariable('rhum', 'f8', ('bulk_time', 'eta_rho', 'xi_rho'))
    var.long_name = 'relative humidity'
    var.units = 'fraction'

    var = nw.createVariable('prate', 'f8', ('bulk_time', 'eta_rho', 'xi_rho'))
    var.long_name = 'precipitation rate'
    var.units = 'cm day-1'

    var = nw.createVariable('wspd', 'f8', ('bulk_time', 'eta_rho', 'xi_rho'))
    var.long_name = 'wind speed 10m'
    var.units = 'm s-1'

    var = nw.createVariable('radlw', 'f8', ('bulk_time', 'eta_rho', 'xi_rho'))
    var.long_name = 'net outgoing longwave radiation'
    var.units = 'Watts meter-2'
    var.positive = 'upward flux, cooling water'

    var = nw.createVariable('radlw_in', 'f8', ('bulk_time', 'eta_rho', 'xi_rho'))
    var.long_name = 'downward longwave radiation'
    var.units = 'Watts meter-2'
    var.positive = 'downward flux, warming water'

    var = nw.createVariable('radsw', 'f8', ('bulk_time', 'eta_rho', 'xi_rho'))
    var.long_name = 'solar shortwave radiation'
    var.units = 'Watts meter-2'
    var.positive = 'downward flux, heating water'

    var = nw.createVariable('sustr', 'f8', ('bulk_time', 'eta_u', 'xi_u'))
    var.long_name = 'surface u-momentum stress'
    var.units = 'Newton meter-2'

    var = nw.createVariable('svstr', 'f8', ('bulk_time', 'eta_v', 'xi_v'))
    var.long_name = 'surface v-momentum stress'
    var.units = 'Newton meter-2'

    var = nw.createVariable('uwnd', 'f8', ('bulk_time', 'eta_u', 'xi_u'))
    var.long_name = '10m u-wind component'
    var.units = 'm/s'

    var = nw.createVariable('vwnd', 'f8', ('bulk_time', 'eta_v', 'xi_v'))
    var.long_name = '10m v-wind component'
    var.units = 'm/s'

    nw.title = title
    nw.date = str(datetime.datetime.now())
    nw.grd_file = grdname
    nw.type = 'CROCO heat flux bulk forcing file'

    for tndx, t in enumerate(bulkt):
        if tndx % 20 == 0:
            print(f'Time Step Bulk: {tndx} of {len(bulkt)}')
        nw.variables['bulk_time'][tndx] = t

    nw.close()
