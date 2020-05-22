import pandas as pd
import numpy as np
import h5py
import os
import matplotlib.pyplot as plt
import argparse
from scipy import optimize
from scipy import stats
import glob
import pygama
from pygama.analysis import histograms
from pygama.analysis import peak_fitting
import json


def main():

    #read tier 2 runs for Ba data
    detector = "I02160A"
    t2_folder = "/lfs/l1/legend/detector_char/enr/hades/char_data/"+detector+"/tier2/ba_HS4_top_dlt/pygama/"
    keys, data = read_all_t2(t2_folder)

    no_events = data.size #all events
    print("No. events: ", no_events)

    key = "e_ftp"
    key_data = obtain_key_data(data, keys, key, no_events)

    no_bins = 10000 #7722=ideal number for 0.5 kev bin width

    #Linearly calibrated data:
    print("")
    print("Linearly calibrating energy...")

    with open('calibration_coef.json') as json_file:
        calibration_coefs = json.load(json_file)
        m = calibration_coefs['m']
        m_err = calibration_coefs['m_err']
        c = calibration_coefs['c']
        c_err = calibration_coefs['c_err']

    print("m: ", m, " , c: ", c)

    calibrated_energy = (key_data-c)/m
    counts, bins_cal, bars = plt.hist(calibrated_energy, bins=no_bins)
    plt.close("all")

    #_________Construct dlt observable________
    print("")
    print("Constructing Ba133 dead layer observable...")



    

def read_all_t2(t2_folder):
    "get data from all tier2 files from same run within a directory"

    run1_files = glob.glob(t2_folder + "/*run0001*.h5")

    file_list = []
    for filename in run1_files:
        df = pd.read_hdf(filename, "data")
        file_list.append(df)
    
    df_total = pd.concat(file_list, axis=0, ignore_index=True)
    keys = df_total.keys()
    data = df_total.to_numpy()

    return keys, data


def obtain_key_data(data, keys, key, no_events):
    "obtain data values for a particular key, e.g. energy"

    key_index = (np.where(keys == key))[0][0]
    key_data = data[:no_events, key_index]

    return key_data

def gaussian(x,a,b,c,d):
    "gaussian function with offset d"
    f = a*np.exp(-((x-b)**2.0/(2.0*c**2.0))) +d
    return f

def gaussian_cdf(x,e,f,g):
    "gaussian cdf function"
    f = g*stats.norm.cdf(x, e, f) #default e=0=mean/loc, f=1=sigma/scale
    return f

def gaussian_and_bkg(x, a, b, c, d, e, f):
    "fit function for 356kev peak"
    f = gaussian(x, a, b, c, d) + gaussian_cdf(x, e, f, g)
    return f


#     gaussian = a*np.exp(-((x-b)**2.0/(2.0*c**2.0))) +d
#     cdf_bkg = 

def linear_fit(x, m, c):
    "linear function"
    f = m*x + c
    return f

def quadratic_fit(x,a,b,c):
    "quadratic function"
    f = a*x**2 + b*x + c
    return f

def sqrt_curve(x,a,c):
    "square root function with offset"
    f = a*np.sqrt(x+c)
    return f

def chi_sq_calc(xdata, ydata, yerr, fittype, popt):
    "calculate chi sq and p-val of a fit given the data points and fit parameters, e.g. fittype ='linear'"
   
    y_obs = ydata
    y_exp = []

    for index, y_i in enumerate(y_obs):
        x_obs = xdata[index]
        if fittype == "linear":
            y_exp_i = linear_fit(x_obs, *popt)
        if fittype == "gaussian":
            y_exp_i = gaussian(x_obs, *popt)
        if fittype == "quadratic":
            y_exp_i = quadratic_fit(x_obs, *popt)
        if fittype == "sqrt":
            y_exp_i = sqrt_curve(x_obs, *popt)
        y_exp.append(y_exp_i)

    #chi_sq, p_value = stats.chisquare(y_obs, y_exp)#this is without errors
    chi_sq = 0.
    residuals = []
    for i in range(len(y_obs)):
        if yerr[i] != 0:
            residual = (y_obs[i]-y_exp[i])/(yerr[i])
        else:
            residual = 0.
        chi_sq += (residual)**2
        residuals.append(residual)

    N = len(y_exp) #number of data points
    dof = N-1
    chi_sq_red = chi_sq/dof

    p_value = 1-stats.chi2.cdf(chi_sq, dof)

    return chi_sq, p_value, residuals, dof


def fit_peak(key, bins, counts, xmin, xmax): #p_guess):
    "fit a gaussian to a peak and return mean and FWHM range"

    no_bins = bins.size 

    xdata = []
    ydata = []
    for bin in bins:
        bin_centre = bin + 0.5*(max(bins)-(min(bins)))/no_bins #this leads to incorrect indexing so leave out for now
        if bin_centre < xmax and bin_centre > xmin:
            xdata.append(bin_centre)
            bin_index = np.where(bins == bin)[0][0]
            ydata.append(counts[bin_index])

    xdata = np.array(xdata)
    ydata = np.array(ydata)     

    yerr = np.sqrt(ydata) #counting error

    #initial rough guess of gaussian params
    aguess = max(ydata) - min(ydata)
    aguess_index = np.where(ydata == max(ydata))[0][0]
    bguess = xdata[aguess_index]
    cguess = (xmax-xmin)/2
    dguess = min(ydata)
    p_guess = [aguess, bguess, cguess, dguess]
    bounds=(0, [np.inf, np.inf, np.inf, np.inf])
    sigma = []
    for index, i in enumerate(yerr):    
        if i != 0:
            sigma.append(yerr[index])
        else:
            sigma.append(1) #just to prevent errors...
    sigma = np.array(sigma)
    popt, pcov = optimize.curve_fit(gaussian, xdata, ydata, p0=p_guess, sigma = sigma, maxfev = 1000000, method ="trf", bounds = bounds)
    mu, sigma = np.abs(popt[1]), np.abs(popt[2]) #must be positive
    mu_err, sigma_err = np.sqrt(pcov[1][1]), np.sqrt(pcov[2][2])
    
    fig, ax = plt.subplots()
    ax.errorbar(xdata, ydata, xerr=0, yerr =yerr, label = "Data", elinewidth = 1, fmt='x', ms = 0.75, mew = 3.0)
    xfit = np.linspace(min(xdata), max(xdata), 1000)
    plt.plot(xfit, gaussian(xfit,*popt), "g", label = "Gaussian fit")
    plt.xlabel(key)
    plt.ylabel("Counts")
    plt.legend()

    chi_sq, p_value, residuals, dof = chi_sq_calc(xdata, ydata, yerr, "gaussian", popt)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    info_str = '\n'.join((r'$\mu=%.2f \pm %.2f$' % (mu, mu_err, ), r'$\sigma=%.2f \pm %.2f$' % (np.abs(sigma), sigma_err,), r'$\chi^2/dof=%.2f/%.0f$'%(chi_sq, dof))) #, r'$p=%.3g$'%p_value))
    ax.text(0.70, 0.20, info_str, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)

    FWHM = 2*np.sqrt(2*np.log(2))*sigma #gaussian FWHM relationship
    peak_range = [mu-FWHM, mu+FWHM]

    return mu, sigma, mu_err, sigma_err


if __name__ =="__main__":
    main()