import numpy as np
import matplotlib.pyplot as plt
from rc_params import plot_params
pp = plot_params()
from astropy import constants
import YY
yy = YY.YonseiYale()
sb = 5.67e-8
L_sun = 3.846e26
R_sun = 6.955e8

def radius_sb(L, T):
    # Stephan-Boltzmann
#     return np.sqrt(L / (constants.sigma_sb * T**4))
    return np.sqrt(L / (sb * T**4))

def luminosity_sb(R, T):
    return constants.sigma_sb * R**2 * T**4

# M, L, T = yy.mass[0], 10**(yy.log_L[0])*L_sun, 10**(yy.log_teff[0])
M, L, T = yy.mass, yy.log_L, 10**(yy.log_teff)
R = radius_sb(L, T) / R_sun

plt.clf()
# plt.plot(M, R)
plt.plot(T, L)
# plt.xlim(plt.gca().get_xlim()[::-1])
plt.xlim(8000, 4000)
plt.ylim(0, 20)
# plt.xlabel("$\mathrm{Mass~(M}_\odot\mathrm{)}$")
# plt.ylabel("$\mathrm{Radius~(R}_\odot\mathrm{)}$")
plt.show()
