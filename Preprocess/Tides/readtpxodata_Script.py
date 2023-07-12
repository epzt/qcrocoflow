import numpy as np
from netCDF4 import Dataset
from matplotlib.path import Path



def readTPXOdata(TPXOfile, ROMSnames, lon, lat):
    """
    This function reads tidal data from the TPXO model and extracts relevant data based on provided ROMS names.

    Parameters:
    TPXOfile (dict): dictionary with TPXO data files (grids and data)
    ROMSnames (list): list of names of ROMS tidal constituents to be extracted
    lon (numpy.array): longitudes of the area of interest
    lat (numpy.array): latitudes of the area of interest

    Returns:
    TPXO (dict): dictionary with extracted data from TPXO model
    """
    # Calculate the min and max bounds for given longitudes and latitudes
    lonmin, lonmax = np.min(lon), np.max(lon)
    latmin, latmax = np.min(lat), np.max(lat)
    # Create bounds for the longitude and latitude
    bndx = [lon[0, 0], lon[-1, 0], lon[-1, -1], lon[0, -1]]
    bndy = [lat[0, 0], lat[-1, 0], lat[-1, -1], lat[0, -1]]
    # Variables to extract from the TPXO file
    variables = ['h', 'U', 'V']
    # Initialize the counter for 'atlas30' grids
    count_atlas30 = -1

    # Initialize the TPXO dictionary to store all data
    TPXO = {var: {'x_a30': None, 'y_a30': None, 'depth_a30': None, 'mask_a30': None, 'z_a30': None} for var in
            variables}
    TPXO['harmonic_a30'] = []

    for harmonic in ROMSnames:
        TPXOfile_grid = TPXOfile['grid']['atlas30']
        count_atlas30 += 1
        TPXO['harmonic_a30'].append({'harmonic': harmonic})

        # Access elevation and velocity data for the current ROMS name
        TPXOfile_elev = TPXOfile['elev'][harmonic.lower()]
        TPXOfile_vel = TPXOfile['vel'][harmonic.lower()]
        files = [TPXOfile_elev, TPXOfile_vel, TPXOfile_vel]

        for i, variable in enumerate(variables):
            coordinate = 'z' if variable == 'h' else variable.lower()
            # Open the data file and extract longitude and latitude data
            data = Dataset(files[i], 'r')
            X1 = data['lon_' + coordinate][:]
            Y1 = data['lat_' + coordinate][:]
            # Create a grid from longitudes and latitudes
            X, Y = np.meshgrid(X1, Y1)

            # Open the TPXO grid file and extract depth data
            with Dataset(TPXOfile_grid, 'r') as nc_file:
                H = nc_file.variables['h' + coordinate][:]
                H = np.transpose(H)

            # Read real and imaginary parts of the variables
            Re_part = np.array(Dataset(files[i], 'r').variables[variable.lower() + 'Re'][:])
            Im_part = np.array(Dataset(files[i], 'r').variables[variable.lower() + 'Im'][:])

            # Create a complex number from the real and imaginary parts
            H_complex = Re_part + 1j * Im_part
            # Transpose the H_complex matrix to match X and Y dimensions
            H_complex = np.transpose(H_complex)

            # Get the indices of cells that are inside the bounds defined earlier
            I, J = np.where((Y >= latmin) & (Y <= latmax) & (X >= lonmin) & (X <= lonmax))

            # Remove duplicate indices
            I, J = np.unique(I), np.unique(J)

            # Create subsets of X, Y and H_complex using found indices
            x = X[I][:, J]
            y = Y[I][:, J]
            z = H_complex[I][:, J] * (0.001 if variable == 'h' else 0.0001)

            # Compute depth
            depth = (H[I][:, J]).T
            depth = np.transpose(depth)

            # Create a mask for points that are inside the bounds defined
            polygon = Path(np.column_stack((bndx, bndy)))
            mask = polygon.contains_points(np.column_stack((x.ravel(), y.ravel()))).reshape(x.shape)

            # Assigning the computed data to the TPXO dictionary
            TPXO[variable]['x_a30'] = x
            TPXO[variable]['y_a30'] = y
            TPXO[variable]['depth_a30'] = depth
            TPXO[variable]['mask_a30'] = mask

            # Initialize the H_complex array in the TPXO dictionary if it is None
            if TPXO[variable]['z_a30'] is None:
                TPXO[variable]['z_a30'] = np.empty((mask.shape[0], mask.shape[1], len(ROMSnames)), dtype=complex)

            # Assign the computed H_complex to the TPXO dictionary

            if variable == 'h':
                TPXO[variable]['z_a30'][:, :, count_atlas30] = z
            else:
                TPXO[variable]['z_a30'][:, :, count_atlas30] = z / np.tile(TPXO[variable]['depth_a30'], (1, 1, 1))



    # Return the dictionary with all the data
    return TPXO
