import numpy as np

def hanning_file(h):
    M, L = h.shape
    Mm = M - 1
    Mmm = M - 2
    Lm = L - 1
    Lmm = L - 2

    h[1:Mm, 1:Lm] = 0.125 * (h[0:Mmm, 1:Lm] + h[2:M, 1:Lm] +
                             h[1:Mm, 0:Lmm] + h[1:Mm, 2:L] +
                             4 * h[1:Mm, 1:Lm])
    h[0, :] = h[1, :]
    h[Mm, :] = h[Mmm, :]
    h[:, 0] = h[:, 1]
    h[:, Lm] = h[:, Lmm]

    return h


