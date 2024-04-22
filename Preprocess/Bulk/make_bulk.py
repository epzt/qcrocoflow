import os
from netCDF4 import Dataset
import numpy as np
import xarray as xr
from .create_nc_bulk import create_bulk
from .nc_add_globatt import add_global_attributes
from .Interp_era5 import interp_ERA5
from PyQt5.QtCore import Qt
from PyQt5.QtCore import Qt
from PyQt5 import QtGui, QtWidgets
def compute_time_values(time, itolap, dt, tlen0):
    """
    Compute the time values for the overlap at the beginning and the end
    """
    # Changer le curseur de la souris à un sablier

    for aa in range(itolap):
        time[aa] = time[itolap] - (itolap - aa) * dt

    for aa in range(itolap):
        time[tlen0 + itolap + aa] = time[tlen0 + itolap] + aa * dt
    return time


def process_files(era5_dir, y_min, y_max, m_min, m_max, grid_dir, title, grdname, Dmin, Yorig, Mmin, directory):
    """
    Process all files in the given directory
    """
    # Changer le curseur de la souris à un sablier
    QtWidgets.QApplication.setOverrideCursor(QtGui.QCursor(Qt.WaitCursor))
    try:
        # Open the croco grid file
        croco_grid = xr.open_dataset(os.path.join(grdname))

        # Extract the variables
        eta_rho = croco_grid['eta_rho']
        xi_rho = croco_grid['xi_rho']
        angle = croco_grid['angle']

        for Y in range(y_min, y_max + 1):
            if Y == y_min:
                mo_min = m_min
            else:
                mo_min = 1

            if Y == y_max:
                mo_max = m_max
            else:
                mo_max = 12

            for M in range(mo_min, mo_max + 1):
                print('\n')
                print(f"Processing year {Y} - month {M}")
                print('\n')

                # Open the NetCDF file
                with Dataset(os.path.join(era5_dir, f"T2M_Y{Y}M{M}.nc")) as nc:

                    # Extract the time variable
                    era5_time = nc.variables['time'][:]

                # Compute the mean difference between consecutive time values
                dt = np.mean(np.gradient(era5_time))
                print(f"dt={dt}")

                # Compute the length of the time variable
                tlen0 = len(era5_time)
                print(f"tlen0={tlen0}")

                # Compute the overlap
                freq = 1  # hourly
                itolap_era5 = 2
                itolap = freq * itolap_era5

                # Compute the total length of the time array
                tlen = tlen0 + 2 * itolap
                print(f"tlen={tlen}")
                print(f"Overlap is {itolap_era5} records before and after")

                # Initialize the time array
                time = np.zeros(tlen)

                # Assign the ERA5 time values to the appropriate portion of the time array
                time[itolap: tlen0 + itolap] = era5_time

                print("=====================")
                print("Compute time for croco file")
                print("=====================")

                time = compute_time_values(time, itolap, dt, tlen0)

                print("=====================")
                print("Create the frc/blk netcdf file")
                print("=====================")
                blkname = f"{title}_bulkY{Y}M{M:02d}.nc"
                directorybulk = os.path.join(directory, blkname)
                print(f"Create a new bulk file: {blkname}")
                create_bulk(grdname, title, time, 0, directorybulk)
                add_global_attributes(directorybulk, Yorig, Mmin, Dmin, 'ERA5')
                interp_ERA5(era5_dir, Y, M, directorybulk, angle, tlen0)
                print(" ")

    finally:
        # Rétablir le curseur de la souris
        QtWidgets.QApplication.restoreOverrideCursor()
def get_y_m(era5_dir):
    """
    Get the minimum and maximum year and month
    """
    # Initialize variables to hold the minimum and maximum year and month
    y_min, y_max, m_min, m_max = float('inf'), float('-inf'), float('inf'), float('-inf')

    # Loop through all the files in the ERA5 directory
    for file in os.listdir(era5_dir):
        # Extract the year and month from the file name
        YM = file.split('_')[1][1:-3]
        Y = int(YM[:4])
        try:
            M = int(YM[4:])
        except ValueError:
            M = int(YM[5:])  # If the month part starts with 'M', skip the 'M'
            print(f"Error converting month to int in file {file}, skipped 'M'")

        # Update the minimum and maximum year and month if necessary
        y_min = min(y_min, Y)
        y_max = max(y_max, Y)
        m_min = min(m_min, M)
        m_max = max(m_max, M)

    return y_min, y_max, m_min, m_max


def main(era5_dir, mo_min, grid_dir, title, grdname, d_min, directory):
    """
    Main function
    """
    # variables
    y_orig = 2000
    mo_min = mo_min
    directory = directory
    d_min = d_min

    grid_dir = grid_dir
    era5_dir = era5_dir
    title = title
    grdname = grdname

    y_min, y_max, m_min, m_max = get_y_m(era5_dir)
    process_files(era5_dir, y_min, y_max, m_min, m_max, grid_dir, title, grdname, d_min, y_orig, mo_min, directory )


if __name__ == "__main__":
    main()
