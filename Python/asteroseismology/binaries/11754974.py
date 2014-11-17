import numpy as np
import matplotlib.pyplot as plt
import pyfits
import glob
from scipy.signal.spectral import lombscargle
from scaling_relations import nu_max, delta_nu
from astero import astero
from rc_params import plot_params
reb, rfb = plot_params()
from colours import plot_colours
ocols = plot_colours()

def testing():
    x = np.arange(0, 10, .1)
    y = np.sin(2*x + np.pi/4.)
    yerr = np.ones_like(y)*.01
    plt.clf()
    plt.plot(x, y)
    ys, A = fit_single_sine_err(x, y, yerr, 2)
    plt.plot(x, ys, 'r')
    plt.show()
    print find_phase(A)/np.pi

# calculate the periodogram for plotting
def pgram(x, y, peak1):
    fs = np.linspace(5, 45, 10000) # c/d
    ws = 2*np.pi*fs  # lombscargle uses angular frequencies
    pgram = lombscargle(x, y, ws)
    plt.clf()
    plt.subplot(2, 1, 1)
#     plt.errorbar(x, y, yerr=yerr, fmt='k.', capsize=0, ecolor='.8')
    plt.plot(x, y, 'k.')
    plt.subplot(2, 1, 2)
    plt.xlabel('$\mu~c/d$')
    plt.axvline(peak1, color='r')
    plt.plot(ws/(2*np.pi), pgram)
    plt.show()
    raw_input('enter')

# find the phases of each segment
def all_phases(xl, yl, yerrl, peak):
    nsegs = len(xl)
    phis = np.zeros(nsegs)
    for i in range(nsegs):
        if len(xl[i]):  # make sure you don't have an empty array
            ysl, Al = fit_single_sine_err(xl[i], yl[i], yerrl[i], peak)
            phis[i] = find_phase(Al)
    return phis

# returns list of lists of the full time series broken down into nsegs
# ndays long
def segment(x, y, yerr, ndays, nsegs):  # nseg is zero counted
    xs, ys, yerrs = [], [], []
    for seg in range(nsegs):
        l = (x[0]+ndays*seg < x) * (x < x[0]+ndays*seg+ndays)
        xs.append(x[l])
        ys.append(y[l])
        yerrs.append(yerr[l])

    BJD = np.zeros(len(xs))
    for i in range(len(xs)):
        if len(xs[i]):
            n = 1
            BJD[i] = xs[i][0]
        else:
            BJD[i] = xs[i-n][0] + 10 * n
            n += 1

    return xs, ys, yerrs, BJD

def light_travel_time(phis, nu):
    phi_mean = np.mean(phis)
    return (phis - phi_mean) / (2*np.pi*nu)

def fit_single_sine_err(x, y, yerr, w):
    M = np.ones((len(x), 2+1))
    S = yerr**2 * np.eye(len(yerr))
    M[:, 0] = np.sin(w*x)
    M[:, 1] = np.cos(w*x)
    A = np.linalg.solve(M.T.dot(np.linalg.inv(S)).dot(M),
                        (M.T.dot(np.linalg.inv(S)).dot(y)))
    ys = np.zeros_like(x)
    ys = A[0]*np.sin(w*x) + A[1]*np.cos(w*x)
    ys += A[2]
    return ys, A

def fit_single_sine(x, y, w):
    M = np.ones((len(x), 2+1))
    M[:, 0] = np.sin(w*x)
    M[:, 1] = np.cos(w*x)
    A = np.linalg.solve(np.dot(M.T, M), np.dot(M.T, y))
    ys = np.zeros_like(x)
    ys = A[0]*np.sin(w*x) + A[1]*np.cos(w*x)
    ys += A[2]
    return ys, A

def find_phase(A):
     return np.arctan(A[1]/A[0])

# load data, median normalise and join together
def load_join(KID, nquarters, sc=False):
    if sc == True:
        lc_files = \
                glob.glob("/Users/angusr/.kplr/data/lightcurves/%s/*_slc.fits"
                          % str(KID).zfill(9))  # short cadence
    else:
        lc_files = \
                glob.glob("/Users/angusr/.kplr/data/lightcurves/%s/*_llc.fits"
                          % str(KID).zfill(9))  # long cadence
    print len(lc_files), 'files found'
    for i, lc_file in enumerate(lc_files):
        hdulist = pyfits.open(lc_file)
        tbdata = hdulist[1].data
        t = np.ascontiguousarray(tbdata["TIME"], dtype=np.float64)
        flux = tbdata["PDCSAP_FLUX"]
        flux_err = tbdata["PDCSAP_FLUX_ERR"]
        q = tbdata["SAP_QUALITY"]
        n = np.isfinite(t)*np.isfinite(flux)*np.isfinite(flux_err)*(q==0)
        t, flux, flux_err = t[n], flux[n], flux_err[n]
        fluxmed = np.median(flux)
        flux /= fluxmed
        flux_err /= fluxmed
        if i == 0:
            x, y, yerr = t, flux, flux_err
        else:
            x = np.concatenate((x, t))
            y = np.concatenate((y, flux))
            yerr = np.concatenate((yerr, flux_err))
        if i == nquarters-1:
            break

    # convert ys to float64
    y2 = np.empty((len(y)))
    for i in range(len(y)):
        y2[i] = y[i].astype('float64')
    y = y2

    return x, y, yerr

if __name__ == "__main__":

    # load data
    KID = "11754974"
    x, y, yerr = load_join(KID, 15)

#     plt.clf()
#     plt.errorbar(x, y, yerr=yerr, **reb)
#     plt.show()
#     raw_input('enter')

    peak1, peak2, peak3 = 16.34, 21.40, 20.91

    # divide data into nsegs segments of ndays days
    ndays = 10
    nsegs = int((max(x)-min(x))/ndays)

    print nsegs, 'nsegs'
    xl, yl, yerrl, BJD = segment(x, y, yerr, ndays, nsegs)

    phis1 = all_phases(xl, yl, yerrl, peak1)
    print phis1, np.mean(phis1)
    times1 = light_travel_time(phis1, peak1) * 24*3600
    phis2 = all_phases(xl, yl, yerrl, peak2)
    times2 = light_travel_time(phis2, peak2) * 24*3600
    phis3 = all_phases(xl, yl, yerrl, peak3)
    times3 = light_travel_time(phis3, peak3) * 24*3600

    times = np.vstack((times1, times2, times3))
    mean_times = np.mean(times, axis=0)

    plt.clf()
#     plt.plot(BJD, phis1, '.')
    plt.plot(BJD, times1, '.')
    plt.plot(BJD, times2, '.')
    plt.plot(BJD, times3, '.')
    plt.plot(BJD, mean_times)
    plt.show()
