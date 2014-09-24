import numpy as np
import matplotlib.pyplot as plt
import pyfits
from scipy.signal.spectral import lombscargle
from scipy.signal import periodogram

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

plt.subplot(2,1,1)
plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')

plt.subplot(2,1,2)

ps = periodogram(y)
plt.plot(ps[0], ps[1], 'b-')
# plt.xlim(0, .1)
plt.ylim(0, 1e7)
plt.xlabel('frequency')

plt.show()
