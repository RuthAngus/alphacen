import numpy as np
import matplotlib.pyplot as plt
import triangle
from params import plot_params
reb = plot_params()
from colours import plot_colours
cols = plot_colours()
import scipy.interpolate as spi

def interp(xold, yold, xnew, s):
    tck = spi.splrep(xold, yold, s=s)
    ynew = spi.splev(xnew, tck, der=0)
    return ynew

f8, f8_err, rho, rho_err = np.genfromtxt("flickers.dat").T
x, xerr = rho, rho_err
y, yerr = f8, f8_err

n1, c1 = np.genfromtxt("CODAchain1.txt").T
n2, c2 = np.genfromtxt("CODAchain2.txt").T

alpha_chain = np.concatenate((c1[1:10000], c2[1:10000]))
alpha = np.median(alpha_chain)
beta_chain = np.concatenate((c1[10001:20000], c2[10001:20000]))
beta = np.median(beta_chain)
sigma_chain = np.concatenate((c1[20001:30000], c2[20001:30000]))
sigma = np.median(sigma_chain)

print alpha, beta, sigma

xm = np.mean(x)
x -= xm

# sample posteriors
nsamp = 10
a = np.random.choice(alpha_chain, nsamp)
b = np.random.choice(beta_chain, nsamp)
s = np.random.choice(sigma_chain, nsamp)

plt.clf()
xs = np.linspace(min(x), max(x), 100)
# l = yerr==max(yerr)
# x, y, xerr, yerr = x[l], y[l], xerr[l], yerr[l]

# for i in range(nsamp):
#     plt.plot(xs, a[i]+b[i]*xs+np.random.randn(1)*(1./s[i]**2)
#              + np.random.randn(1)*yerr, color=".5", alpha=.1)
#     for j in range(len(x)):
#         plt.plot(np.ones(nsamp)*x[j],
#                  a[i]+b[i]*x[j]+np.random.randn(nsamp)*(1./s[i]**2)
#                  + np.random.randn(nsamp)*yerr[j], color=cols.blue, alpha=.05)

# for i in range(nsamp):
#     plt.plot(xs, a[i]+b[i]*xs+np.random.randn(1)*(1./s[i]**2),
#              color=".5", alpha=.1)

inds = np.argsort(x)
x = x[inds]
y = y[inds]
yerr = yerr[inds]

ys = alpha + beta * x
tau = 1./sigma**2
xs = np.linspace(min(x), max(x), 1000)
tau_tot = tau + yerr

# intrinsic scatter
plt.fill_between(xs+xm, alpha+beta*xs-tau, alpha+beta*xs+tau, color=cols.blue,
        alpha=.5, edgecolor="w")

# intrinsic scatter + individual uncertainties, 1 sigma
ys1 = interp(x, ys-tau_tot, xs, 1)
ys2 = interp(x, ys+tau_tot, xs, 1)
plt.fill_between(xs+xm, ys1, ys2, color=cols.blue, alpha=.2, edgecolor="w")

# intrinsic scatter + individual uncertainties, 2 sigma
ys1 = interp(x, ys-2*tau_tot, xs, 4.5)
ys2 = interp(x, ys+2*tau_tot, xs, 4.5)
plt.fill_between(xs+xm, ys1, ys2, color=cols.blue, alpha=.2, edgecolor="w")

plt.errorbar(x+xm, y, xerr=xerr, yerr=yerr, fmt="k.", capsize=0, alpha=.5,
             ecolor=".5", mec=".2")
plt.plot(xs+xm, alpha+beta*xs, ".2", linewidth=1)

plt.xlim(min(xs)+xm, max(xs)+xm)
plt.xlabel("$\log_{10}(\\rho_{star}[\mathrm{kgm}^{-3}])$")
plt.ylabel("$\log_{10}\mathrm{(8-hour~flicker~[ppm])}$")
plt.savefig("flicker")

# samples = np.vstack((alpha_chain, beta_chain, sigma_chain)).T
# fig = triangle.corner(samples)
# fig.savefig("triangle")
