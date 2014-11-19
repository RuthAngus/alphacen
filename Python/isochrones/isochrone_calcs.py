import numpy as np
import matplotlib.pyplot as plt
from rc_params import plot_params
pp = plot_params()
from astropy import constants
from YREC import YREC
yrec = YREC()

L_sun = 3.846e26
R_sun = 6.955e8

def radius_sb(L, T):
    L *= L_sun
    sb = 5.67e-8
    return np.sqrt(L / (sb * T**4)) / R_sun

# for a given mass and teff, what is the radius?
def rad(M, T):

    # load isochrones
    iso_M, iso_L, iso_T = yrec.mass, 10**yrec.logL, yrec.teff

    # find closest mass, closest teff and corresponding l
    m = min(iso_M, key=lambda x:abs(x-M))
    teff = min(iso_T, key=lambda x:abs(x-T))
    l = iso_L[iso_M==m]  # warning: this is often degenerate!

    return radius_sb(l[0], teff), m, teff, l

if __name__ == "__main__":

    rad(1.38375, 6000.)

#     plt.clf()
#     plt.plot(T, yrec.logL)
#     plt.xlim(8000, 4000)
