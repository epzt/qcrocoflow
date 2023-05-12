import numpy as np

def get_angle_file(latu, lonu, spheroid="wgs84"):
    if spheroid[:3] == "sph":
        A = 6371000.0
        B = A
        E = np.sqrt(A*A - B*B) / A
        EPS = E*E / (1 - E*E)
    elif spheroid[:3] == "cla":
        A = 6378206.4
        B = 6356583.8
        E = np.sqrt(A*A - B*B) / A
        EPS = E*E / (1 - E*E)
    elif spheroid[:3] == "iau":
        A = 6378160
        B = 6356774.516
        E = np.sqrt(A*A - B*B) / A
        EPS = E*E / (1 - E*E)
    elif spheroid[:3] == "wgs":
        A = 6378137
        E = 0.081819191
        B = np.sqrt(A**2 - (A*E)**2)
        EPS = E*E / (1 - E*E)
    else:
        raise ValueError("Unknown spheroid specified!")

    latu = np.radians(latu)
    lonu = np.radians(lonu)

    latu[latu == 0] = np.finfo(float).eps

    M, L = latu.shape

    PHI1 = latu[:, :-1]
    XLAM1 = lonu[:, :-1]
    PHI2 = latu[:, 1:]
    XLAM2 = lonu[:, 1:]

    PHI2[PHI1 == PHI2] = PHI2[PHI1 == PHI2] + 1e-14
    XLAM2[XLAM1 == XLAM2] = XLAM2[XLAM1 == XLAM2] + 1e-14

    xnu1 = A / np.sqrt(1.0 - (E * np.sin(PHI1))**2)
    xnu2 = A / np.sqrt(1.0 - (E * np.sin(PHI2))**2)

    TPSI2 = (1 - E*E) * np.tan(PHI2) + E*E * xnu1 * np.sin(PHI1) / (xnu2 * np.cos(PHI2))

    DLAM = XLAM2 - XLAM1
    CTA12 = (np.cos(PHI1) * TPSI2 - np.sin(PHI1) * np.cos(DLAM)) / np.sin(DLAM)
    azim = np.arctan(1 / CTA12)

    DLAM2 = (np.abs(DLAM) < np.pi) * DLAM + (DLAM >= np.pi) * (-2 * np.pi + DLAM) + (DLAM <= -np.pi) * (2 * np.pi + DLAM)
    azim = azim + (azim < -np.pi) * 2 * np.pi - (azim >= np.pi) * 2 * np.pi
    azim = azim + np.pi * np.sign(-azim) * (np.sign(azim) != np.sign(DLAM2))
    angle = np.zeros((M, L+1))
    angle[:, 1:L] = (np.pi / 2) - azim
    angle[:, 0] = angle[:, 1]
    angle[:, -1] = angle[:, -2]

    return angle




