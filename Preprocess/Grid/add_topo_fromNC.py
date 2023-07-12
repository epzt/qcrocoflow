import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata

def add_topo_file(grdname, topofile):
    """
    This function interpolates the topography data into the CROCO grid.

    Parameters:
    grdname: str
        Name of the CROCO grid file.
    topofile: str
        Name of the topography file.

    Returns:
    h: ndarray
        Interpolated topography data.
    """

    # Load CROCO grid
    try:
        grd = Dataset(grdname, 'r')
    except FileNotFoundError:
        print(f"File {grdname} not found.")
        return None

    lon_rho = grd.variables['lon_rho'][:]
    lat_rho = grd.variables['lat_rho'][:]
    pm = grd.variables['pm'][:]
    pn = grd.variables['pn'][:]
    grd.close()

    # Compute the average resolution of the CROCO grid
    dx = np.mean(1 / pm)
    dy = np.mean(1 / pn)
    dx_croco = np.mean([dx, dy])

    # Define the bounding box for subsetting TOPO data
    dl = max([1, 2 * (dx_croco / (60 * 1852))]) # estimate the required padding
    bbox = [np.min(lon_rho) - dl, np.max(lon_rho) + dl, np.min(lat_rho) - dl, np.max(lat_rho) + dl]

    # Load TOPO file
    try:
        nc = Dataset(topofile, 'r')
    except FileNotFoundError:
        print(f"File {topofile} not found.")
        return None

    lon_topo = nc.variables['lon'][:]
    lat_topo = nc.variables['lat'][:]

    # Subset the TOPO data based on the bounding box
    # The longitude is adjusted to the [-180, 180] range
    lon_topo_subset = lon_topo[np.logical_and(lon_topo >= bbox[0], lon_topo <= bbox[1])]
    lat_topo_subset = lat_topo[np.logical_and(lat_topo >= bbox[2], lat_topo <= bbox[3])]

    # Retrieve the topo data corresponding to the subset
    topo_subset = -nc.variables['topo'][np.ix_(lat_topo_subset, lon_topo_subset)]
    nc.close()

    # Rescale the resolution of the topo subset data
    dx_topo = compute_avg_resolution(lon_topo_subset, lat_topo_subset)
    while dx_croco > dx_topo:
        lon_topo_subset, lat_topo_subset, topo_subset = downscale_resolution(lon_topo_subset, lat_topo_subset, topo_subset)
        dx_topo = compute_avg_resolution(lon_topo_subset, lat_topo_subset)

    # Interpolate the topo data into the CROCO grid using cubic interpolation
    h = griddata((lon_topo_subset.ravel(), lat_topo_subset.ravel()), topo_subset.ravel(), (lon_rho, lat_rho), method='cubic')

    return h

def compute_avg_resolution(lon, lat):
    """
    Compute the average resolution of the grid defined by lon and lat.

    Parameters:
    lon: ndarray
        Longitudes.
    lat: ndarray
        Latitudes.

    Returns:
    dx_avg: float
        Average resolution.
    """
    R = 6367442.76 # Earth radius in meters
    deg2rad = np.pi / 180
    dx = R * deg2rad * np.diff(lon) * np.cos(deg2rad * lat[:-1])
    dy = R * deg2rad * np.diff(lat)
    dx_avg = np.mean(np.hypot(dx, dy))

    return dx_avg

def downscale_resolution(lon, lat, data):
    """
    Downscale the resolution of the data by a factor of 2.

    Parameters:
    lon: ndarray
        Longitudes.
    lat: ndarray
        Latitudes.
    data: ndarray
        Data.

    Returns:
    lon: ndarray
        Downscaled longitudes.
    lat: ndarray
        Downscaled latitudes.
    data: ndarray
        Downscaled data.
    """
    lon = 0.5 * (lon[1:] + lon[:-1])
    lon = lon[::2]
    lat = 0.5 * (lat[1:] + lat[:-1])
    lat = lat[::2]
    data = 0.25 * (data[1:, :-1] + data[1:, 1:] + data[:-1, :-1] + data[:-1, 1:])
    data = data[::2, ::2]

    return lon, lat, data
