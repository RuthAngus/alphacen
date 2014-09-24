import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const

# My model starts with a time of transit (t90), period, and mass.
# I solve for time of periapse:

# the true anomaly, f is defined as 90 or pi/2 at the time of central transit.
# the time dependence comes in because you know the period and you know the time of central transit.

Mstar = .770 * const.M_sun
Mstar_err = 0.018 * const.M_sun
ecc = 0.
w = 0.
f = 0. # no idea what this should be for a circular orbit!!!
M2 = .1 * const.M_sun
a = .1 * const.au

def rv_model(f, Mstar, Mstar_err, ecc, w, M2, a):

    f90 = np.pi/2-w
    E90 = np.arctan2(np.sqrt(1-ecc**2)*np.sin(f90),(ecc+np.cos(f90)))
    M90 = E90 - ecc*np.sin(E90)

    return np.sqrt(const.G/((Mstar + M2)*a*const.au*(1-ecc**2)))\
            * M2*(np.cos(w+f)+ecc*np.cos(w))/100.0

#     return (t90 - M90/2/np.pi*period) # find time of periapsis

print rv_model(f, Mstar, Mstar_err, ecc, w, M2, a)

# Mstar = 0.770 \pm 0.018 M_sun, a is derived from Kepler's third law,
# f is the true anomaly, ecc is still zero. M2 is the planet mass,
# and G and AU are constants to make the units work out.

# def lnlike(x):
#
#     v_t = rv_model(f, Mstar, Mstar_err, ecc, w)
#     return np.sum(np.log((1-x[-1])/(err_t**2 + x[-4]**2)**0.5 *
#                   np.exp(-0.5*((v_t + x[-5] - rv_t)**2/(err_t**2 + x[-4]**2)))
#                   + x[-1]/(err_t**2 + x[-4]**2 + x[-3]**2)**0.5 *
#                   np.exp(-0.5*(v_t + x[-5] - rv_t)**2/(err_t**2 + x[-4]**2 + x[-3]**2))))
#
# x = Pb, S, Sp, V0

# before adding the prior.

# x[-1] is the probabilitiy a point is bad, x[-4] is our jitter, x[-3] is our excess jitter. x[-5] is the RV zeropoint.

# rv_t is our observations, v_t is our (circular) model.

# We are assuming a transit center of 2454833.0 + 1865.13 \pm 0.01 days.
