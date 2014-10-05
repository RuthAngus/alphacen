import numpy as np
import matplotlib.pyplot as plt
import pyfits
from scipy.signal.spectral import lombscargle
from scipy.signal import periodogram
import george
from george.kernels import ExpSquaredKernel, ExpSine2Kernel, WhiteKernel, CosineKernel
from scipy.misc import derivative
from scipy.optimize import minimize, fmin

def predict(xs, theta, y, yerr):
    k = theta[0] * ExpSquaredKernel(theta[1]) * ExpSine2Kernel(theta[2], theta[4])
    k += WhiteKernel(theta[3])
    gp = george.GP(k)
    print type(yerr), np.shape(yerr), type(y), np.shape(y)
    gp.compute(xs, yerr)
    return gp.predict(y, xs)[0]

def fit_rv(theta, x, y, yerr, xs):

    # check for normalisation
    if np.median(y) > .5:
        print "not normalised!"

    k = theta[0] * ExpSquaredKernel(theta[1]) * ExpSine2Kernel(theta[2], theta[4])
    k += WhiteKernel(theta[3])
    gp = george.GP(k)
    gp.compute(x, yerr)
    print "initial likelihood = ", gp.lnlikelihood(y, quiet=True)

    mu, cov = gp.predict(y, xs)

    # plot initial guess
    plt.clf()
    plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')
    plt.plot(xs, mu, 'r')
    plt.savefig('init_rv_fit')

#     pars, res = gp.optimize(x, y, yerr)
#     results = fmin(neglnlike, theta, args=(x, y, yerr))
#     print results
#     print pars, 'pars'
#     raw_input('enter')

#     k = pars[0] * ExpSquaredKernel(pars[1])  # * ExpSine2Kernel(pars[2], pars[4])
#     k += WhiteKernel(pars[2])
#     gp = george.GP(k)
#     gp.compute(x, yerr)
#     print "final likelihood = ", gp.lnlikelihood(y, quiet=True)

#     mu, cov = gp.predict(y, xs)

    derivs = derivative(predict, xs, args=(theta, y, yerr))

    return mu, cov, derivs
