import numpy as np

def get_angle_file(latu, lonu, spheroid="wgs84"):
    """Calculate the angle between geodetic and spherical coordinates.

    Args:
        latu (np.array): 2D array of latitudes.
        lonu (np.array): 2D array of longitudes.
        spheroid (str, optional): The name of the reference ellipsoid to use.
            Default is "wgs84".

    Returns:
        angle (np.array): 2D array of calculated angles.

    """

    # Constants for different spheroids
    SPHEROIDS = {
        'sph': (6371000.0, 6371000.0),
        'cla': (6378206.4, 6356583.8),
        'iau': (6378160, 6356774.516),
        'wgs': (6378137, 6378137 * (1 - 0.081819191**2)**0.5)
    }

    # Get constants for the given spheroid
    A, B = SPHEROIDS.get(spheroid[:3])
    if not A or not B:
        raise ValueError("Unknown spheroid specified!")
    E = np.sqrt(A*A - B*B) / A


    # Convert degrees to radians
    latu = np.radians(latu)
    lonu = np.radians(lonu)

    # Replace zero values to prevent division by zero
    latu[latu == 0] = np.finfo(float).eps

    # Get size of the input arrays
    M, L = latu.shape

    # Split the arrays into two parts
    PHI1, XLAM1 = latu[:, :-1], lonu[:, :-1]
    PHI2, XLAM2 = latu[:, 1:], lonu[:, 1:]

    # Add small values to identical coordinates to prevent division by zero
    PHI2[PHI1 == PHI2] += 1e-14
    XLAM2[XLAM1 == XLAM2] += 1e-14

    # Calculate intermediate variables
    xnu1 = A / np.sqrt(1.0 - (E * np.sin(PHI1))**2)
    xnu2 = A / np.sqrt(1.0 - (E * np.sin(PHI2))**2)
    TPSI2 = (1 - E*E) * np.tan(PHI2) + E*E * xnu1 * np.sin(PHI1) / (xnu2 * np.cos(PHI2))
    DLAM = XLAM2 - XLAM1
    CTA12 = (np.cos(PHI1) * TPSI2 - np.sin(PHI1) * np.cos(DLAM)) / np.sin(DLAM)

    # Calculate azimuth and correct its range
    azim = np.arctan(1 / CTA12)
    DLAM2 = (np.abs(DLAM) < np.pi) * DLAM + (DLAM >= np.pi) * (-2 * np.pi + DLAM) + (DLAM <= -np.pi) * (2 * np.pi + DLAM)
    azim += (azim < -np.pi) * 2 * np.pi - (azim >= np.pi) * 2 * np.pi
    azim += np.pi * np.sign(-azim) * (np.sign(azim) != np.sign(DLAM2))

    # Calculate the final angles
    angle = np.zeros((M, L+1))
    angle[:, 1:L] = (np.pi / 2) - azim
    angle[:, 0] = angle[:, 1]
    angle[:, -1] = angle[:, -2]

    return angle
