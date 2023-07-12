from scipy.interpolate import griddata
import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator


def interpTPXO(TPXOvar, count, lon, lat, mask):
    """
    This function extracts TPXO data and interpolates it onto CROCO grid.
    Inputs:
        - TPXOvar: hRE, hIm, uRE, uIm, vRe, or vIm
        - count:
        - lon, lat: 2D arrays of longitude and latitude to interpolate to
        - mask: array of 1s and 0s corresponding to wet and dry points of lon & lat
    Output:
        - VARinterp (VAR on lon, lat grid)
    """

    # Initialize the interpolated variable as a complex array of zeros
    VARinterp = np.zeros(mask.shape, dtype=np.complex128)
    iswet = mask == 1  # Array of booleans where the mask is 1 (wet points)

    # Extract coordinates and variable values from TPXO data
    x = TPXOvar['x_a30']
    y = TPXOvar['y_a30']
    z = np.squeeze(TPXOvar['z_a30'][:, :, count])
    depth_masked = TPXOvar['depth_a30']
    depth = depth_masked.filled()  # Convert the MaskedArray to an ndarray, filling masked values with a default fill value

    m = np.logical_and(TPXOvar['mask_a30'], depth > 0)  # Logical mask where depth is above 0 and TPXOvar mask is true

    # Interpolate z using nearest griddata method
    z = griddata((x[m], y[m]), z[m], (x, y), method='nearest')

    # TODO: debug
    points = np.array([x, y])

    # Création de la grille d'interpolation
    F = griddata((x.flatten(), y.flatten()), z.flatten(), (lon[iswet].T, lat[iswet].T), method='linear')

    # Assignation des valeurs interpolées aux positions correspondantes dans VARinterp
    VARinterp[iswet] = F
    VARinterp = VARinterp.transpose()

    return VARinterp
