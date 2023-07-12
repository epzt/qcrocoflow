def rho2uvp(rfield):
    """
    This function converts the rho coordinates to the corresponding u, v and p.

    Parameters:
        rfield (np.array): The rho-coordinates to be converted.

    Returns:
        ufield (np.array): The corresponding u-coordinates.
        vfield (np.array): The corresponding v-coordinates.
        pfield (np.array): The corresponding p-coordinates.
    """

    # Get the shape of the rho field
    Mp, Lp = rfield.shape

    # Calculate the size of the u, v, and p grids
    M = Mp - 1
    L = Lp - 1

    # Convert the rho field to the v field
    vfield = 0.5 * (rfield[0:M, :] + rfield[1:Mp, :])

    # Convert the rho field to the u field
    ufield = 0.5 * (rfield[:, 0:L] + rfield[:, 1:Lp])

    # Convert the u field to the p field
    pfield = 0.5 * (ufield[0:M, :] + ufield[1:Mp, :])

    return ufield, vfield, pfield
