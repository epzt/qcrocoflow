import numpy as np


def ap2ep(Au, PHIu, Av, PHIv):
    """
    This function converts amplitude and phase tidal constituents into ellipse parameters.

    Parameters:
    Au (numpy.array): amplitude of u-component
    PHIu (numpy.array): phase of u-component
    Av (numpy.array): amplitude of v-component
    PHIv (numpy.array): phase of v-component

    Returns:
    SEMA (numpy.array): Semi-major axis of the ellipse
    ECC (numpy.array): ECC of the ellipse
    INC (numpy.array): INC of the ellipse
    PHA (numpy.array): Phase angle of the ellipse
    """
    # Convert the phases from degrees to radians
    PHIu = np.deg2rad(PHIu)
    PHIv = np.deg2rad(PHIv)

    # Calculate the complex amplitudes for u and v
    i = complex(0, 1)
    u = Au * np.exp(-i * PHIu)
    v = Av * np.exp(-i * PHIv)

    # Compute the counterclockwise and clockwise circles
    wp = (u + i * v) / 2
    wm = np.conj(u - i * v) / 2

    # Calculate the amplitudes and angles
    Wp = np.abs(wp)
    Wm = np.abs(wm)
    THETAp = np.angle(wp)
    THETAm = np.angle(wm)

    # Compute the ellipse parameters
    SEMA = Wp + Wm  # Semi-major axis
    SEMI = Wp - Wm  # Semi-minor axis
    ECC = np.where(SEMA != 0, SEMI / SEMA, 0)  # ECC
    PHA = (THETAm - THETAp) / 2  # Phase angle
    INC = (THETAm + THETAp) / 2  # INC

    # Convert back from radians to degrees
    PHA = np.rad2deg(PHA)
    INC = np.rad2deg(INC)

    # Adjust phase and INC angles to range [0, 2*pi)
    PHA[PHA < 0] += 360
    INC[INC < 0] += 360

    print('hello')
    return SEMA, ECC, INC, PHA
