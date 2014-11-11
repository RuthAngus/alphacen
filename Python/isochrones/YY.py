import numpy as np
import matplotlib.pyplot as plt
import glob

class YonseiYale(object):

    def __init__(self):

        x, z = '74', '01'
        x, z = '59', '06'

        DIR = '/Users/angusr/angusr/isochrones/a0o2'
        files = glob.glob('%s/x%sz%s/m*x%sz%s.track1' % (DIR, x, z, x, z))

        npts = 24
        ages = np.zeros((npts, len(files)))
        log_teffs = np.zeros((npts, len(files)))
        log_Ls = np.zeros((npts, len(files)))
        masses = np.zeros((npts, len(files)))

        ms = []
        for i, f in enumerate(files):
            m = np.float('%s.%s' % (f[45:46], f[46:47]))

            # load isochrones
            ages[:, i], log_teffs[:, i], log_Ls[:, i] = \
                    np.genfromtxt(f, skip_header=1, usecols=(1, 2, 3)).T
            masses[:, i] = np.ones(npts)*m

        self.age = ages
        self.log_teff = log_teffs
        self.log_L = log_Ls
        self.mass = masses

if __name__ == "__main__":

    YY = YonseiYale()
    print YY.age[0]
    print YY.log_teff[0]
    print YY.log_L[0]
    print YY.mass[0]
