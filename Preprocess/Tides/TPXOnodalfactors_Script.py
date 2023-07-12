import numpy as np
from datetime import datetime


def compute_factors(t):
    """
    function to compute different factors for f and u.
    """

    VN = np.mod(360 * (0.719954 - 5.372617 * t + 0.000006 * t * t), 360)
    if VN < 0:
        VN += 360
    VN = VN * np.pi / 180

    cN = np.cos(VN)
    c2N = np.cos(2 * VN)
    c3N = np.cos(3 * VN)
    sN = np.sin(VN)
    s2N = np.sin(2 * VN)
    s3N = np.sin(3 * VN)

    return cN, c2N, c3N, sN, s2N, s3N


def TPXOnodalfactors(dnum, names):
    """
    Computes f and u factors for the harmonics listed in cell array names.
    Only those which are in the TPXO model are evaluated.
    They are evaluated at time dnum.
    """

    f = dict()
    u = dict()
    f_factors = dict()
    u_factors = dict()
    days_until_1900 = ((datetime(1900, 1, 1) - datetime(1, 1, 1)).days + 367)

    t = ((dnum + 0.5) - days_until_1900) / 36525

    cN, c2N, c3N, sN, s2N, s3N = compute_factors(t)

    f["mm"] = 1.0 - 0.1300 * cN + 0.0013 * c2N
    u["mm"] = 0
    f["mf"] = 1.0429 + 0.4135 * cN - 0.004 * c2N
    u["mf"] = -0.4143 * sN + 0.0468 * s2N - 0.0066 * s3N
    f["o1"] = 1.0089 + .1871 * cN - 0.0147 * c2N + 0.0014 * c3N
    u["o1"] = 0.1885 * sN - 0.0234 * s2N + 0.0033 * s3N
    f["k1"] = 1.0060 + 0.1150 * cN - 0.0088 * c2N + 0.0006 * c3N
    u["k1"] = -0.1546 * sN + 0.0119 * s2N - 0.0012 * s3N
    f["m2"] = 1.0004 - 0.0373 * cN + 0.0002 * c2N
    u["m2"] = -0.0374 * sN
    f["k2"] = 1.0241 + 0.2863 * cN + 0.0083 * c2N - 0.0015 * c3N
    u["k2"] = -0.3096 * sN + 0.0119 * s2N - 0.0007 * s3N

    f["q1"] = f["o1"]
    u["q1"] = u["o1"]
    f["n2"] = f["m2"]
    u["n2"] = u["m2"]
    f["mn4"] = f["m2"] ** 2
    u["mn4"] = 2 * u["m2"]
    f["m4"] = f["m2"] ** 2
    u["m4"] = 2 * u["m2"]
    f["ms4"] = f["m2"]
    u["ms4"] = u["m2"]

    # Assign f and u factors for each harmonic name
    for n in names:
        if n in f:
            f_factors[n] = f[n]
            u_factors[n] = np.mod(u[n] * 180.0 / np.pi, 360.0)
        else:
            f_factors[n] = 1
            u_factors[n] = 0

    return f_factors, u_factors