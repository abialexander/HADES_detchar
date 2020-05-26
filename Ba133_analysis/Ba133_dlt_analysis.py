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
    print("Available keys: " ,keys)



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

    xmin_356, xmax_356 = 353, 360 #kev
    plt.figure()
    popt, pcov, xfit = fit_peak_356("Energy (keV)", bins_cal, counts, xmin_356, xmax_356)
    counts, bins, bars = plt.hist(calibrated_energy, bins=no_bins, histtype='step', color='grey')
    plt.xlim(xmin_356, xmax_356) 
    plt.ylim(100, 10**6)
    # plt.yscale('log')
    # plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/plots/356keV_dlt.png")

    # plt.figure()
    plt.plot(xfit, -1*popt[1]*gaussian_cdf(xfit,popt[6],popt[7]) + popt[2], "r--", label ="-b*gauss_cdf(x,g,h) + c (bkg)")
    plt.yscale("log")
    plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/plots/356keV_dlt.png")



    

def read_all_t2(t2_folder):
    "get data from all tier2 files from same run within a directory"

    run1_files = glob.glob(t2_folder + "/*run0001*.h5")

    file_list = []
    for filename in run1_files:
        df = pd.read_hdf(filename, "data")
        file_list.append(df)
    
    df_total = pd.concat(file_list, axis=0, ignore_index=True)

    #print(df_total['data']['e_ftp']) #this doesnt work

    keys = df_total.keys()
    data = df_total.to_numpy()

    return keys, data


def obtain_key_data(data, keys, key, no_events):
    "obtain data values for a particular key, e.g. energy"

    key_index = (np.where(keys == key))[0][0]
    key_data = data[:no_events, key_index]

    return key_data

def gaussian(x,a,b,c):
    "gaussian function without offset"
    f = a*np.exp(-((x-b)**2.0/(2.0*c**2.0)))
    return f

def gaussian_cdf(x,a,b):
    "gaussian cdf function"
    f = stats.norm.cdf(x, a, b) #default e=0=mean/loc, f=1=sigma/scale
    return f

def gaussian_and_bkg(x, a, b, c, d, e, f, g, h):
    "fit function for 356kev peak"
    f = a*gaussian(x, d, e, f) - b*gaussian_cdf(x, g, h) + c
    return f

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
        if fittype == "gaussian_and_bkg":
            y_exp_i = gaussian_and_bkg(x_obs, *popt)
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

def fit_peak_356(key, bins, counts, xmin, xmax):
    "fit the 356 keV peak with gaussian +cdf bkg"

    no_bins = bins.size 

    xdata = []
    ydata = []
    for bin in bins:
        bin_centre = bin + 0.5*(max(bins)-(min(bins)))/no_bins 
        if bin_centre < xmax and bin_centre > xmin:
            xdata.append(bin_centre)
            bin_index = np.where(bins == bin)[0][0]
            ydata.append(counts[bin_index])

    xdata = np.array(xdata)
    ydata = np.array(ydata)     

    yerr = np.sqrt(ydata) #counting error

    #initial rough guess of params
    aguess = max(ydata) - min(ydata) #gauss amplitude
    aguess_index = np.where(ydata == max(ydata))[0][0]
    bguess = min(ydata)/100 #cdf amp
    cguess = 0 #min(ydata) #offset
    dguess = 1 #intrinsic gauss amp
    eguess = xdata[aguess_index] #gauss mean
    fguess = 1 #gauss sigma
    gguess = eguess #cdf mean
    hguess = 1 #fguess #cdf sigma
    # bguess = xdata[aguess_index] #gauss mean
    # cguess = (xmax-xmin)/2 #gauss sigma
    # dguess = min(ydata) #offset
    # eguess = bguess #cdf mean
    # fguess = cguess #cdf sigma
    # gguess = aguess/100 #cdf amplitude
    p_guess = [aguess, bguess, cguess, dguess, eguess, fguess, gguess, hguess]
    #bounds=([0, 0, -np.inf, 0, 0, 0, 0, 0], [np.inf]*8)
    sigma = []
    for index, i in enumerate(yerr):    
        if i != 0:
            sigma.append(yerr[index])
        else:
            sigma.append(1) #just to prevent errors...
    sigma = np.array(sigma)
    popt, pcov = optimize.curve_fit(gaussian_and_bkg, xdata, ydata, p0=p_guess, sigma = sigma, maxfev = 1000000, method ="trf") #, bounds = bounds)
    
    a,b,c,d,e,f,g,h = popt[0],popt[1],[2],popt[3],popt[4],popt[5],popt[6],popt[7] 

    fig, ax = plt.subplots()
    #ax.errorbar(xdata, ydata, xerr=0, yerr =yerr, label = "Data", elinewidth = 1, fmt='x', ms = 0.75, mew = 3.0)
    plt.errorbar(xdata, ydata, xerr=0, yerr =yerr, label = "Data", elinewidth = 1, fmt='x', ms = 0.75, mew = 3.0)
    xfit = np.linspace(min(xdata), max(xdata), 1000)
    plt.plot(xfit, gaussian_and_bkg(xfit,*popt), "g", label = "a*gauss(x,d,e,f) - b*gauss_cdf(x,g,h) + c") 
    plt.plot(xfit, -1*b*gaussian_cdf(xfit,g,h) + c, "r--", label ="-b*gauss_cdf(x,g,h) + c")
    plt.xlabel(key)
    plt.ylabel("Counts")
    plt.legend(loc="upper right", fontsize=8)

    chi_sq, p_value, residuals, dof = chi_sq_calc(xdata, ydata, yerr, "gaussian_and_bkg", popt)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    info_str = '\n'.join((r'$a=%.3g \pm %.3g$' % (popt[0], np.sqrt(pcov[0][0])), r'$b=%.3g \pm %.3g$' % (popt[1], np.sqrt(pcov[1][1])), r'$c=%.3g \pm %.3g$' % (popt[2], np.sqrt(pcov[2][2])), r'$d=%.3g \pm %.3g$' % (popt[3], np.sqrt(pcov[3][3])), r'$e=%.3g \pm %.3g$' % (popt[4], np.sqrt(pcov[4][4])), r'$f=%.3g \pm %.3g$' % (popt[5], np.sqrt(pcov[5][5])),r'$g=%.3g \pm %.3g$' % (popt[6], np.sqrt(pcov[6][6])), r'$h=%.3g \pm %.3g$' % (popt[7], np.sqrt(pcov[7][7])), r'$\chi^2/dof=%.2f/%.0f$'%(chi_sq, dof)))
    plt.text(0.02, 0.98, info_str, transform=ax.transAxes, fontsize=8,verticalalignment='top', bbox=props) #ax.text..ax.tra

    return popt, pcov, xfit


if __name__ =="__main__":
    main()