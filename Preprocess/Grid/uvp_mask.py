import numpy as np

def uvp_mask_file(rfield):
    """
    This function computes the u, v, and p fields from the rho mask by taking the
    product of adjacent elements in the horizontal and vertical directions.

    Parameters:
        rfield (np.array): The input r-field array.

    Returns:
        ufield (np.array): The computed u-field array.
        vfield (np.array): The computed v-field array.
        pfield (np.array): The computed p-field array.
    """

    # Determine the shape of rfield
    Mp, Lp = rfield.shape

    # Compute the u-field by multiplying horizontally adjacent elements in rfield
    ufield = np.zeros((Mp, Lp-1))
    ufield[:, :] = rfield[:, :-1] * rfield[:, 1:]

    # Compute the v-field by multiplying vertically adjacent elements in rfield
    vfield = np.zeros((Mp-1, Lp))
    vfield[:, :] = rfield[:-1, :] * rfield[1:, :]

    # Compute the p-field by multiplying vertically and horizontally adjacent elements in rfield
    pfield = np.zeros((Mp-1, Lp-1))
    pfield[:, :] = rfield[:-1, :-1] * rfield[:-1, 1:] * rfield[1:, :-1] * rfield[1:, 1:]

    return ufield, vfield, pfield
