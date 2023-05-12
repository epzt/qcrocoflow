import numpy as np

def spheric_dist_file(lat1, lat2, lon1, lon2):
    # Earth radius
    R = 6367442.76

    # Determine proper longitudinal shift.
    l = np.abs(lon2 - lon1)
    l[l >= 180] = 360 - l[l >= 180]

    # Convert decimal degrees to radians.
    deg2rad = np.pi / 180
    lat1 = lat1 * deg2rad
    lat2 = lat2 * deg2rad
    l = l * deg2rad

    # Compute the distances
    dist = R * np.arcsin(np.sqrt(((np.sin(l) * np.cos(lat2))**2) + (((np.sin(lat2) * np.cos(lat1)) - (np.sin(lat1) * np.cos(lat2) * np.cos(l)))**2)))
    return dist
