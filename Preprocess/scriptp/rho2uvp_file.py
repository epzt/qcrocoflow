def rho2uvp(rfield)
    Mp, Lp = rfield.shape
    M = Mp - 1
    L = Lp - 1

    vfield = 0.5  (rfield[0M, ] + rfield[1Mp, ])
    ufield = 0.5  (rfield[, 0L] + rfield[, 1Lp])
    pfield = 0.5  (ufield[0M, ] + ufield[1Mp, ])

    return ufield, vfield, pfield
