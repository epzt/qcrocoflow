import numpy as np


def qsat(Ta, Pa=None):
    """Computes specific humidity at saturation.

    Args:
        Ta: Air temperature in degrees Celsius.
        Pa: Air pressure in millibars (optional).

    Returns:
        q: Saturation specific humidity in kg/kg.
    """
    Pa = 1015

    # Calculate saturation vapor pressure in millibars
    ew = 6.1121 * (1.0007 + 3.46e-6 * Pa) * np.exp((17.502 * Ta) / (240.97 + Ta))

    # Calculate specific humidity at saturation
    q = 0.62197 * (ew / (Pa - 0.378 * ew))

    return q
