import numpy as np
import matplotlib.pyplot as plt
import glob

class YREC(object):

    def __init__(self):

        DIR = '/Users/angusr/angusr/isochrones/YREC'

        ages = [1.122, 1.259, 1.413, 1.585, 1.778, 1.995, 2.239, 2.512, 2.818,
                3.162, 3.548, 3.981, 4.467, 5.012, 5.623, 6.310, 7.079, 7.943,
                8.913, 10.0, 11.220, 12.589, 14.125, 15.849]

        mass, teff, logL, logg = \
                np.genfromtxt('%s/p_000_corr.txt' % DIR,
                              skip_header=2, usecols=(0, 1, 2, 3)).T

        self.mass = mass
        self.teff = teff
        self.logL = logL
        self.logg = logg

if __name__ == "__main__":

    yrec = YREC()

    plt.clf()
    plt.plot(yrec.teff, 10**yrec.logL, 'k.')
    plt.ylim(0, 10)
    plt.xlim(8000, 4000)
    plt.show()
