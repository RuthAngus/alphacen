import numpy as np
import matplotlib.pyplot as plt
import orbit
import emcee
import triangle
import nestle
import george
from george.kernels import ExpSquaredKernel, ExpSine2Kernel

# Simple one planet model:
def oneplanetmodel(theta, t):
    P1, K1, T01, V01, Ecc1, omega1 = theta
    return orbit.radvel(t, P1, K1, T01, V01, Ecc1, omega1)

# The likelihood function:
def lnlike(theta, t, rv_obs, rverr):
    model_rv = oneplanetmodel(theta, t)
    chisq = np.sum(((rv_obs - model_rv) / rverr)**2)
    return -chisq / 2.

# restrictive priors
def lnprior(theta):
    if 100. < theta[0] < 300. and 100 < theta[1] < 300 \
            and 2455200. < theta[2] < 2455300. and 0. < theta[3] < 10e6 \
            and 0. < theta[4] < .1 and 0. < theta[5] < 90.:
                return 0.
#     if 0. < theta[0] < 500. and 0. < theta[1] < 500 \
#             and 2455000. < theta[2] < 2455500. and 0. < theta[3] < 10e6 \
#             and 0. < theta[4] < .2 and 0. < theta[5] < 90.:
#                 return 0.
    return -np.inf

def lnprob(theta, t, rv_obs, rverr):
    return lnprior(theta) + lnlike(theta, t, rv_obs, rverr)

# FECH/CHIRON data
t1, rv1, rv_err1 = np.genfromtxt('/Users/angusr/Python/rvexp/FCHIP105854.txt',
                              skip_header=1).T

plt.clf()
plt.errorbar(t1, rv1, yerr=rv_err1, fmt='k.', capsize=0, ecolor='.8')

# load FEROZ data
t2, rv2, rv_err2 = np.genfromtxt('/Users/angusr/Python/rvexp/FEROZHIP105854.txt',
                              skip_header=1).T

plt.errorbar(t2, rv2, yerr=rv_err2, fmt='b.', capsize=0, ecolor='.8')

t = np.concatenate((t1, t2))
rv = np.concatenate((rv1, rv2))
rv_err = np.concatenate((rv_err1, rv_err2))

p, p_err, K, K_err, a, a_err, e, e_err, w, w_err, m, m_err, T0, t0_err, g1, \
        g1_err, g2, g2_err = np.genfromtxt('HIP105854_params.txt', skip_header=1).T

V0 = 0.
theta = p, K, T0, V0, e, w

x = np.linspace(min(t), max(t), 1000)
model_rvs = oneplanetmodel(theta, x)

plt.plot(x, model_rvs)
plt.savefig('initial')

print 'burn in...'
nwalkers, ndim = 32, len(theta)
p0 = [theta+1e-4*np.random.rand(ndim) for i in range(nwalkers)]
args = [t, rv, rv_err]
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=args)
p0, lp, state = sampler.run_mcmc(p0, 2000)

print 'production run...'
sampler.reset()
p0, lp, state = sampler.run_mcmc(p0, 4000)

fig_labels = ['P', 'K', 'T0', 'V0', 'Ecc', 'omega']
flatchain = sampler.chain[:, 50:, :].reshape((-1, ndim))
fig = triangle.corner(flatchain, truths=theta, labels=fig_labels)
plt.savefig('triangleHIP')

mcmc_result = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                  zip(*np.percentile(flatchain, [16, 50, 84], axis=0)))
mres = np.array(mcmc_result)[:, 0]

plt.clf()
plt.errorbar(t, rv, yerr=rv_err, fmt='k.', capsize=0, ecolor='.8')
plt.plot(x, oneplanetmodel(mres, x))
plt.savefig('final')
