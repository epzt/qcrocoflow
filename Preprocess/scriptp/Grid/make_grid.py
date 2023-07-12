"""
This module provides a function to create a regular grid from the coordinates
minimum and maximum latitude and longitude provided.
"""

# Python built-in modules
import os
# Third-party libraries
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc_module
from netCDF4 import Dataset
from qgis.core import (QgsCoordinateReferenceSystem, QgsFeature, QgsGeometry,
                       QgsMapLayer, QgsPointXY, QgsProject, QgsRectangle,
                       QgsVectorLayer, QgsRasterLayer, QgsCoordinateTransform)
# Local application/library specific imports
from .bathy_view import BathyView
from .add_topo_from_raster import add_topo_file2
from .create_grid import create_grid_file
from .add_topo_fromNC import add_topo_file
from ...qcrocoflow_qcrocotools import netCDFtoRaster
from .get_rho_script import get_rho
from .uvp_mask import uvp_mask_file
from .get_angle import get_angle_file
from .rho2uvp_script import rho2uvp


def make_grid_function(title, grdname, lon_min, lon_max, lat_min, lat_max, dl, topofile, directory):
    """
        Creates a regular grid from latitude and longitude coordinates
        minimum and maximum provided.

        Parameters:
        title (str): The title of the output file.
        grdname (str): The name of the output file.
        lon_min (float): The minimum longitude of the area of interest.
        lon_max (float): The maximum longitude of the area of interest.
        lat_min (float): The minimum latitude of the area of interest.
        lat_max (float): The maximum latitude of the area of interest.
        dl (float): The interval of latitude and longitude coordinates.
        topofile (str): The path to the topographic file to use for the interpolation.
        directory (str): The directory where the output file will be saved.

        Returns:
        grdname_full_path
        """

    # Create longitude and latitude grids with dl resolution.
    x = np.arange(lon_min, lon_max, dl)
    y = np.arange(lat_min, lat_max, dl)

    # Initialize the coordinate transformation from EPSG:3857 to EPSG:4326.
    transform_4326 = QgsCoordinateTransform(QgsCoordinateReferenceSystem('EPSG:3857'),
                                            QgsCoordinateReferenceSystem('EPSG:4326'),
                                            QgsProject.instance())

    # Create 2D longitude and latitude grids.
    Lonr, Latr = np.meshgrid(x, y)

    # Transform each grid point to EPSG:4326.
    for i in range(Lonr.shape[0]):
        for j in range(Lonr.shape[1]):
            point_4326 = transform_4326.transform(QgsPointXY(Lonr[i, j], Latr[i, j]))
            Lonr[i, j] = point_4326.x()
            Latr[i, j] = point_4326.y()

    # Initialize new variables.
    dndx = np.zeros(Latr.shape)
    dmde = np.zeros(Latr.shape)
    pn = np.full_like(Latr, 1/dl)
    pm = np.full_like(Latr, 1/dl)

    # Convert longitude and latitude grids from rho to u, v and p.
    Lonu, Lonv, Lonp = rho2uvp(Lonr)
    Latu, Latv, Latp = rho2uvp(Latr)

    # Get grid dimensions.
    M, L = Latp.shape

    # Create output NC with specified name and directory.
    grdname_full_path = create_grid_file(L, M, grdname, title, directory)

    # Open netCDF file for writing.
    with nc_module.Dataset(grdname_full_path, 'a') as nc:
        # Assign computed grids to corresponding variables in netCDF file.
            nc['lat_u'][:] = Latu
            nc['lon_u'][:] = Lonu
            nc['lat_v'][:] = Latv
            nc['lon_v'][:] = Lonv
            nc['lat_rho'][:] = Latr
            nc['lon_rho'][:] = Lonr
            nc['lat_psi'][:] = Latp
            nc['lon_psi'][:] = Lonp

    # Initialize matrices for x and y coordinates on the rho grid.
    xr,yr = get_rho(dl)

    # Convert coordinate grids from rho to u, v and p.
    xu, xv, xp = rho2uvp(xr)
    yu, yv, yp = rho2uvp(yr)

    # Calculate the angle.
    angle = get_angle_file(Latu, Lonu)

    # Compute the Coriolis parameter.
    f = 4 * np.pi * np.sin(np.pi * Latr / 180) * 366.25 / (24 * 3600 * 365.25)

    # Fill the grid file with the calculated parameters.
    print('Fill the grid file...')
    with Dataset(grdname_full_path, 'a') as nc:
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

        ######################################
        # Adding topography to the grid      #
        ######################################
        # Identifying the file extension of the topography file
        _, file_extension = os.path.splitext(topofile)
        file_extension = file_extension.lower()
        print(f"The file extension is: {file_extension}")

        # Choosing the topography adding method based on the file extension
        if file_extension == '.nc':
            h = add_topo_file(grdname_full_path, topofile)
        elif file_extension in ['.grd', '.tif']:
            h = add_topo_file2(topofile, x.shape, y.shape)
        else:
            raise ValueError(f"Unsupported file type: {file_extension}")

        #####################################################################
        # Opening the graphical interface to visualize the bathymetry       #
        #####################################################################
        bathy_view = BathyView()
        bathy_view.set_h(h, nc)
        bathy_view.exec_()

        # Retrieve the new depth level
        new_h = bathy_view.save_new_level()

        # Defining the mask based on the depth
        maskr = new_h > 0 if new_h is not None else h > 0

        # Visualizing the mask before processing
        mappable = plt.imshow(maskr, origin='lower', cmap='binary')
        plt.title('Mask before processing')
        plt.colorbar(mappable)
        plt.show()

        # Converting the rho mask to u, v, and p
        masku, maskv, maskp = uvp_mask_file(maskr)

        # Adding the masks to the netCDF file
        with nc_module.Dataset(grdname_full_path, 'a') as nc:
            nc.variables['mask_u'][:] = masku
            nc.variables['mask_v'][:] = maskv
            nc.variables['mask_psi'][:] = maskp
            nc.variables['mask_rho'][:] = maskr

            # Extracting latitude, longitude and "h" data
            lon = nc.variables['lon_rho'][:]
            lat = nc.variables['lat_rho'][:]
            h = nc.variables['h'][:]

        # Name of the GeoTIFF file to create
        tif_filename = "grid.tif"

        # Create an instance of netCDFtoRaster
        nc2r = netCDFtoRaster(QgsCoordinateReferenceSystem(4326))  # Use EPSG:4326 for CRS

        # Create a GeoTIFF file for "h" variable
        nc2r.createRaster(lon, lat, h, tif_filename, h.ndim)

        # Add the GeoTIFF file to the QGIS project
        raster_layer = QgsRasterLayer(tif_filename, "grid_show")

        if not raster_layer.isValid():
            print("Failed to load the layer!")
        else:
            # Add the layer to your QGIS project
            QgsProject.instance().addMapLayer(raster_layer)

        return grdname_full_path