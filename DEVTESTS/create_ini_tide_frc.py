import geopandas as gpd
from scipy.interpolate import RegularGridInterpolator

import numpy as np
import xarray as xr

ncfile = '/home/epoizot/recherche/PROGRAMMES/RIN_seline/Stockline2/Gnbb2/LADEROUTE.nc'
ds = xr.open_dataset(ncfile)
h = ds['h'].to_dataframe().reset_index()

# Load your GeoDataFrame
gdf = gpd.GeoDataFrame(h, geometry=gpd.points_from_xy(h.lon, h.lat), crs='EPSG:4326')

# Extract coordinates and values
coordinates = gdf[['x', 'y']].values
values = gdf['h'].values
print(values)
'''
# Create a grid of points
x = np.linspace(coordinates[:, 0].min(), coordinates[:, 0].max(), 100)
y = np.linspace(coordinates[:, 1].min(), coordinates[:, 1].max(), 100)
X, Y = np.meshgrid(x, y)

# Perform interpolation
Z = RegularGridInterpolator((x, y), values)(np.array([X, Y]))
'''
