import numpy as np
import matplotlib.pyplot as plt
import pyfits
from scipy.signal.spectral import lombscargle
from scipy.signal import periodogram
import george
from george.kernels import ExpSquaredKernel, ExpSine2Kernel, WhiteKernel, CosineKernel
from rc_params import plot_params
from colors import plot_colors
ocols = plot_colors()

lc_file = "/Users/angusr/.kplr/data/lightcurves/005955122/kplr005955122-2009166044711_slc.fits"
hdulist = pyfits.open(lc_file)
tbdata = hdulist[1].data
x = tbdata["TIME"]
y = tbdata["PDCSAP_FLUX"]
yerr = tbdata["PDCSAP_FLUX_ERR"]
q = tbdata["SAP_QUALITY"]
x -= x[0]
n = np.isfinite(x)*np.isfinite(y)*np.isfinite(yerr)*(q==0) * (x < 1)
x, y, yerr = x[n], y[n], yerr[n]
y -= np.median(y)

# from BGdata import BetaGem
# BG = BetaGem()
# x = BG.fHJD - BG.fHJD[0]
# y, yerr = BG.flux-np.median(BG.flux), BG.flux_err
# l = x < .05
# x, y, yerr = x[l], y[l], yerr[l]
# x = np.arange(len(y))

plt.clf()
plt.subplot(2,1,1)
plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')

plt.subplot(2,1,2)
t = x
sp = np.fft.fft(y)
freq = np.fft.fftfreq(t.shape[-1])
print freq
plt.plot(freq, sp.real)
plt.xlim(0, .6)
# plt.xlim(0, 10)
plt.ylim(0, 60000)
# plt.ylim(0, 30000)
plt.show()
raw_input('enter')

theta = [2e5, .05, 10, 1000, .1]
# theta = [.16, .1 ** 2, .5, .01, .5]
xs = np.linspace(min(x), max(x), 100)

# k = kernel(theta)
k = theta[0] * ExpSquaredKernel(theta[1]) * ExpSine2Kernel(theta[2], theta[4])
k += WhiteKernel(theta[3])
gp = george.GP(k)
gp.compute(x, np.sqrt(theta[3]+yerr**2))
print gp.lnlikelihood(y, quiet=True)

print theta
mu, cov = gp.predict(y, xs)
pars, results = gp.optimize(x, y, yerr)
print pars

k = theta[0] * ExpSquaredKernel(theta[1]) * ExpSine2Kernel(theta[2], theta[4])
k += WhiteKernel(theta[3])
gp = george.GP(k)
gp.compute(x, np.sqrt(theta[3]+yerr**2))

nmu, ncov = gp.predict(y, xs)

plt.plot(xs, mu, 'b-')
plt.plot(xs, nmu, color=ocols.blue)
plt.ylabel('$\mathrm{Flux}$')
plt.show()
raw_input('enter')

def predict(xs):
    return gp.predict(y, xs)[0]

from scipy.misc import derivative
derivs = derivative(predict, xs)
print predict(xs)
print derivs

plt.subplot(2,1,2)
rvx = xs
rv_err = .5 + np.random.rand(len(rvx))
# rv_err = .00000001
rvs = derivs/100. #+ np.random.uniform(-1, 1)*rv_err
rvs += np.random.randn(len(rvx)) * rv_err
plt.errorbar(xs, rvs, yerr=rv_err, fmt='k.',
             capsize=0, ecolor='.8')
plt.plot(rvx, derivs/100., color=ocols.pink)
plt.xlabel('$\mathrm{Time~(Days)}$')
plt.ylabel('$\mathrm{RV~(ms}_{-1}\mathrm{)}$')
plt.savefig('rvs')
