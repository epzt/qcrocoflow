import numpy as np

# Importation de la fonction hanning du fichier hanning.py
from .hanning import hanning_file

def smoothgrid_file(h, maskr, hmin, hmax, n_filter_deep_topo, n_filter_final, hmax_coast, r_max):
    def uvp_mask(maskr):
        masku = np.zeros(maskr.shape)
        maskv = np.zeros(maskr.shape)
        maskp = np.zeros(maskr.shape)

        masku[:, 1:] = maskr[:, 1:] * maskr[:, :-1]
        maskv[1:, :] = maskr[1:, :] * maskr[:-1, :]
        maskp[1:, 1:] = maskr[1:, 1:] * maskr[:-1, 1:] * maskr[1:, :-1] * maskr[:-1, :-1]

        return masku, maskv, maskp

    def log_topo_filter(h, maskr, masku, maskv, maskr_ext, hmin, hmax_coast, r_max):
        OneEights = 1/8
        OneThirtyTwo = 1/32

        r = rfact(h, masku, maskv)
        cff = 1.4
        r_max = r_max / cff
        i = 0

        while r > (r_max * cff):
            i += 1

            Lgh = np.log(h / hmin)
            Lgh[hmin == 0] = 0

            lgr_max = np.log((1. + r_max) / (1. - r_max))
            lgr1_max = lgr_max * np.sqrt(2)

            grad = (Lgh[:, 1:] - Lgh[:, :-1])
            cr = np.abs(grad)
            FX = grad * (1. - lgr_max / cr)
            FX[cr <= lgr_max] = 0

            grad = (Lgh[1:, 1:] - Lgh[:-1, :-1])
            cr = np.abs(grad)
            FX1 = grad * (1. - lgr1_max / cr)
            FX1[cr <= lgr1_max] = 0

            grad = (Lgh[1:, :] - Lgh[:-1, :])
            cr = np.abs(grad)
            FE = grad * (1. - lgr_max / cr)
            FE[cr <= lgr_max] = 0

            grad = (Lgh[1:, :-1] - Lgh[:-1, 1:])
            cr = np.abs(grad)
            FE1 = grad * (1. - lgr1_max / cr)
            FE1[cr <= lgr1_max] = 0

            Lgh[1:-1, 1:-1] = Lgh[1:-1, 1:-1] + \
                OneEights * (FX[1:-1, 1:] - FX[1:-1, :-1] + FE[1:, 1:-1] - FE[:-1, 1:-1]) + \
                OneThirtyTwo * (FX1[1:, 1:] - FX1[:-1, :-1] + FE1[1:, :-1] - FE1[:-1, 1:])

            Lgh[0, :] = Lgh[1, :]
            Lgh[-1, :] = Lgh[-2, :]
            Lgh[:, 0] = Lgh[:, 1]
            Lgh[:, -1] = Lgh[:, -2]

            h = hmin * np.exp(Lgh)

            h[(maskr_ext < 0.5) & (h > hmax_coast)] = hmax_coast

            r = rfact(h, masku, maskv)

            if i % 20 == 0:
                print(f"   {i} iterations - r_max = {r}")

        print(f"   {i} iterations - r_max = {r}")

        return h

    def rfact(h, masku, maskv):
        # Calculate gradients
        gradx = h[:, 1:] - h[:, :-1]
        grady = h[1:, :] - h[:-1, :]

        # Calculate ratio of gradient to bathymetry
        rx = np.abs(gradx) / (h[:, 1:] + h[:, :-1])
        ry = np.abs(grady) / (h[1:, :] + h[:-1, :])

        # Apply mask
        rx = rx * masku
        ry = ry * maskv

        # Find the maximum ratio value
        r = max(np.max(rx), np.max(ry))

        return r


    def hanning_smoother(h):
        M, L = h.shape
        Mm = M - 1
        Mmm = M - 2
        Lm = L - 1
        Lmm = L - 2

        h_smoothed = np.copy(h)
        h_smoothed[1:Mm, 1:Lm] = 0.125 * (h[0:Mmm, 1:Lm] + h[2:M, 1:Lm] +
                                         h[1:Mm, 0:Lmm] + h[1:Mm, 2:L] +
                                         4 * h[1:Mm, 1:Lm])

        h_smoothed[0, :] = h_smoothed[1, :]
        h_smoothed[Mm, :] = h_smoothed[Mmm, :]
        h_smoothed[:, 0] = h_smoothed[:, 1]
        h_smoothed[:, Lm] = h_smoothed[:, Lmm]

        return h_smoothed

    def hanning_smoother_coef2d(h, coef):
        # Obtenir les dimensions de la matrice h
        M, L = h.shape
        Mm = M - 1
        Mmm = M - 2
        Lm = L - 1
        Lmm = L - 2
        
        # Initialiser la matrice h_smoothed en copiant h
        h_smoothed = np.copy(h)
        
        # Appliquer le lissage avec les coefficients sur les éléments internes de la matrice
        h_smoothed[1:Mm, 1:Lm] = coef[1:Mm, 1:Lm] * (h[0:Mmm, 1:Lm] + h[2:M, 1:Lm] +
                                                     h[1:Mm, 0:Lmm] + h[1:Mm, 2:L]) \
                                + (1 - 4 * coef[1:Mm, 1:Lm]) * h[1:Mm, 1:Lm]

        # Mettre à jour les bords de la matrice avec les valeurs des voisins
        h_smoothed[0, :] = h_smoothed[1, :]
        h_smoothed[Mm, :] = h_smoothed[Mmm, :]
        h_smoothed[:, 0] = h_smoothed[:, 1]
        h_smoothed[:, Lm] = h_smoothed[:, Lmm]

        return h_smoothed

    def hann_window(h):
        # Définition des constantes de pondération pour le filtre Hanning
        OneFours = 1/4
        OneEights = 1/8
        OneSixteens = 1/16

        # Récupération des dimensions de la matrice h
        M, L = h.shape
        Mm = M - 1
        Mmm = M - 2
        Lm = L - 1
        Lmm = L - 2

        # Création d'une copie de la matrice h pour stocker les valeurs lissées
        h_smoothed = np.copy(h)

        # Application du filtre Hanning 2D sur les éléments internes de la matrice h
        h_smoothed[1:Mm, 1:Lm] = (OneFours * h_smoothed[1:Mm, 1:Lm] +
                                  OneEights * (h_smoothed[0:Mmm, 1:Lm] + h_smoothed[2:M, 1:Lm] +
                                               h_smoothed[1:Mm, 0:Lmm] + h_smoothed[1:Mm, 2:L]) +
                                  OneSixteens * (h_smoothed[0:Mmm, 0:Lmm] + h_smoothed[2:M, 2:L] +
                                                 h_smoothed[0:Mmm, 2:L] + h_smoothed[2:M, 0:Lmm]))

        # Mise à jour des bords de la matrice h_smoothed pour correspondre à ceux de la matrice h
        h_smoothed[0, :] = h_smoothed[1, :]
        h_smoothed[Mm, :] = h_smoothed[Mmm, :]
        h_smoothed[:, 0] = h_smoothed[:, 1]
        h_smoothed[:, Lm] = h_smoothed[:, Lmm]

        # Retourner la matrice h lissée
        return h_smoothed

    print(" Filter topography (new scheme) ...")

    masku, maskv, maskp = uvp_mask(maskr)
    maskr_ext = hann_window(maskr)
    maskr_ext[maskr_ext < 1] = 0

    h[h < hmin] = hmin
    h[h > hmax] = hmax

    if n_filter_deep_topo >= 1:
        print("Apply a filter on the Deep Ocean to reduce isolated seamounts:")
        print(f"   {n_filter_deep_topo} pass of a selective filter.")

        # Construire un coefficient de lissage qui est une fonction linéaire
        # d'une bathymétrie lisse.
        coef = h.copy()
        for i in range(8):
            coef = hann_window(coef)  # coef est une bathymétrie lissée

        coef = 0.125 * (coef / np.max(coef))  # redimensionner la bathymétrie lissée

        for i in range(n_filter_deep_topo):
            h = hanning_smoother_coef2d(h, coef)  # lisser avec coef disponible
            h[(maskr_ext < 0.5) & (h > hmax_coast)] = hmax_coast


    h = log_topo_filter(h, maskr, masku, maskv, maskr_ext, hmin, hmax_coast, r_max)

    if n_filter_final > 1:
    #Ce code s'exécutera si n_filter_final est supérieur à 1. Il lisse la topographie une dernière fois pour prévenir le bruit 2DX. 
    #Ensuite, il met à jour la variable h en appliquant la condition hmin.
        print("Smooth the topography a last time to prevent 2DX noise:")
        print(f"   {n_filter_final} pass of a hanning smoother.")
        
        for i in range(n_filter_final):
            h = hann_window(h)
            h[(maskr_ext < 0.5) & (h > hmax_coast)] = hmax_coast

    # Mettre à jour la variable h en appliquant la condition hmin
    h[h < hmin] = hmin


    return h
