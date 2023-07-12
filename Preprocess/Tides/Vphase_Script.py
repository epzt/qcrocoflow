import numpy as np
from datetime import datetime


def compute_phase(t, DN_Nbr, DN_Pha, DN_Spd, tHour):
    """
    Helper function to compute phase factors.
    """

    Vs = np.mod(360 * (0.751206 + 1336.855231 * t - 0.000003 * t ** 2), 360)
    Vh = np.mod(360 * (0.776935 + 100.002136 * t + 0.000001 * t * t), 360)
    Vp = np.mod(360 * (0.928693 + 11.302872 * t - 0.000029 * t * t), 360)
    VN = np.mod(360 * (0.719954 - 5.372617 * t + 0.000006 * t * t), 360)
    Vp1 = np.mod(360 * (0.781169 + 0.004775 * t + 0.000001 * t * t), 360)

    for V in [Vs, Vh, Vp, VN, Vp1]:
        if V < 0:
            V += 360

    Vdeg = tHour * DN_Spd + np.dot(Vs, DN_Nbr[1]) + np.dot(Vh, DN_Nbr[2]) + \
           np.dot(Vp, DN_Nbr[3]) + np.dot(VN, DN_Nbr[4]) + np.dot(Vp1, DN_Nbr[5]) + DN_Pha

    Vdeg = np.mod(Vdeg, 360)

    return Vdeg


def Vphase(dnum, DN_List):
    """
    Computes equilibrium phases in accordance with Cartwright "tidal analysis - a retrospect".
    """

    days_until_1900 = (datetime(1900, 1, 1) - datetime(1, 1, 1)).days + 367

    t = (dnum + 0.5 - days_until_1900) / 36525

    tHour = np.mod(dnum, 1) * 24
    DN_Nbr = DN_List[:7]
    DN_Pha = DN_List[6]
    DN_Spd = DN_List[7]

    Vdeg = compute_phase(t, DN_Nbr, DN_Pha, DN_Spd, tHour)

    print(Vdeg)
    return Vdeg
