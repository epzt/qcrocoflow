import numpy as np

def process_mask_file(maskin):
    maskout = maskin.copy()
    M, L = maskout.shape
    Mm, Lm = M - 1, L - 1
    Mmm, Lmm = Mm - 1, Lm - 1
    
    neibmask = np.zeros_like(maskout)
    neibmask[1:Mm, 1:Lm] = maskout[:Mmm, 1:Lm] + maskout[2:, 1:Lm] + maskout[1:Mm, :Lmm] + maskout[1:Mm, 2:]
    
    water_neighbor_threshold = 3
    land_neighbor_threshold = 1

    while np.sum(((neibmask[1:Mm, 1:Lm] >= water_neighbor_threshold) & (maskout[1:Mm, 1:Lm] == 0)) | ((neibmask[1:Mm, 1:Lm] <= land_neighbor_threshold) & (maskout[1:Mm, 1:Lm] == 1))) > 0:
        maskout[(neibmask >= water_neighbor_threshold) & (maskout == 0)] = 1
        maskout[(neibmask <= land_neighbor_threshold) & (maskout == 1)] = 0

        maskout[0, 1:Lm] = maskout[1, 1:Lm]
        maskout[Mm, 1:Lm] = maskout[Mm - 1, 1:Lm]
        maskout[1:Mm, 0] = maskout[1:Mm, 1]
        maskout[1:Mm, Lm] = maskout[1:Mm, Lm - 1]

        maskout[0, 0] = min(maskout[0, 1], maskout[1, 0])
        maskout[Mm, 0] = min(maskout[Mm, 1], maskout[Mm - 1, 0])
        maskout[0, Lm] = min(maskout[0, Lm - 1], maskout[1, Lm])
        maskout[Mm, Lm] = min(maskout[Mm, Lm - 1], maskout[Mm - 1, Lm])

        neibmask[1:Mm, 1:Lm] = maskout[:Mmm, 1:Lm] + maskout[2:, 1:Lm] + maskout[1:Mm, :Lmm] + maskout[1:Mm, 2:]

    maskout[:, 0] = maskout[:, 1]
    maskout[:, Lm] = maskout[:, Lm - 1]
    maskout[0, :] = maskout[1, :]
    maskout[Mm, :] = maskout[Mm - 1, :]

    return maskout

