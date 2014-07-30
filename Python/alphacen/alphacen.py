import numpy as np
import pyfits
import matplotlib.pyplot as pl
from mylnlikefn import predict

data = np.genfromtxt("data.txt", skip_header=2).T

t = data[0]
t -= t[0]
rv = data[1]
rv_err2 = data[2]
rv_err = data[3]
bis = data[5]
fwhm = data[6]
logrhk = data[7]
logrhk_err = data[8]

# truncate

# subtract linear trend
ply = np.polyfit(t, rv, 1)
p = np.polyval(ply, t)
rv -= p

# convert rv to m/s
rv *= 1000
rv_err *=1000

# RV = vr*Gp + vc*G
# logrhk = lc*G
# bis = br*Gp + bc*G

# GP regression on the logrhk
xp = np.linspace(min(t), max(t), 100)

# median normalise
logrhk = logrhk/np.median(logrhk) - 1
logrhk_err = 2*logrhk_err/(max(logrhk)-min(logrhk))
logrhk = 2*logrhk/(max(logrhk)-min(logrhk))
logrhk = logrhk-np.median(logrhk)

theta = [0, 0, 0, 0, 0]
# theta = [-2., -2., -1.2, 1., 1.]

yp = predict(xp, t, logrhk, logrhk_err, theta)

# test
t = np.linspace(0, 15, 10)
logrhk = np.sin(t)
logrhk_err = 0.01
xp = np.linspace(min(t), max(t), 100)
yp = predict(xp, t, logrhk, logrhk_err, theta)
print yp

pl.clf()
pl.errorbar(t, logrhk, yerr=logrhk_err, fmt='k.')
# pl.plot(xp, yp, 'b-')
pl.savefig('logrhk')
