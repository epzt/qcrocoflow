import numpy as np

def uvp_mask_file(rfield):
    Mp, Lp = rfield.shape
    M = Mp - 1
    L = Lp - 1
    
    # Calculer vfield en multipliant les éléments adjacents verticalement dans rfield
    vfield = rfield[:-1, :] * rfield[1:, :]
    
    # Calculer ufield en multipliant les éléments adjacents horizontalement dans rfield
    ufield = rfield[:, :-1] * rfield[:, 1:]
    
    # Calculer pfield en multipliant les éléments adjacents verticalement dans ufield
    pfield = ufield[:-1, :] * ufield[1:, :]

    return vfield, ufield, pfield
