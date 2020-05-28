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
    print("356 peak:")


    #fit peak with gaussian and cdf
    xmin_356, xmax_356 = 353, 360 #kev
    plt.figure()
    popt, pcov, xfit = fit_peak_356("Energy (keV)", bins_cal, counts, xmin_356, xmax_356)
    a,b,c,d,e,f,g = popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],popt[6] 
    counts, bins, bars = plt.hist(calibrated_energy, bins=no_bins, histtype='step', color='grey')
    plt.xlim(xmin_356, xmax_356) 
    plt.ylim(100, 10**6)
    plt.plot(xfit, -1*d*gaussian_cdf(xfit,e,f) + g, "r--", label ="-d*gauss_cdf(x,e,f) + g")
    plt.yscale("log")
    plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/plots/356keV_dlt.png")


    #plot signal only and calculate integral
    fig, ax = plt.subplots()
    plt.plot(xfit, gaussian(xfit,a,b,c), "b--", label ="gauss(x,a,b,c)")
    plt.yscale("log")
    integral = gauss_count(a, b, c)
    print("gauss count: ", integral)
    plt.xlim(xmin_356, xmax_356) 
    plt.xlabel("Energy (keV)")
    plt.ylabel("Counts")
    plt.legend(loc="upper right", fontsize=9)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    info_str = r'$C_{356} = %.3g$' % (integral)
    plt.text(0.02, 0.07, info_str, transform=ax.transAxes, fontsize=9,verticalalignment='top', bbox=props) #ax.text
    plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/plots/356keV_dlt_signalonly.png")


    #try other fits - constrained cdf - doesnt work?
    plt.figure()
    popt, pcov, xfit = fit_peak_356_2("Energy (keV)", bins_cal, counts, xmin_356, xmax_356)
    a,b,c,d,e = popt[0],popt[1],[2],popt[3],popt[4]
    counts, bins, bars = plt.hist(calibrated_energy, bins=no_bins, histtype='step', color='grey')
    plt.xlim(xmin_356, xmax_356) 
    plt.ylim(100, 10**6)
    plt.plot(xfit, -1*d*gaussian_cdf(xfit,b,c) + e, "r--", label ="-d*gauss_cdf(x,b,c) + e")
    plt.yscale("log")
    plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/plots/356keV_dlt_2.png")

    

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
    #f = a*np.exp(-(pow((x-b),2)/(2.0*pow(c,2))))
    return f

def gaussian_cdf(x,a,b):
    "gaussian cdf function"
    f = stats.norm.cdf(x, a, b) #default e=0=mean/loc, f=1=sigma/scale
    return f

def gaussian_and_bkg(x, a, b, c, d, e, f, g):
    "fit function for 356kev peak"
    f = gaussian(x, a, b, c) - d*gaussian_cdf(x, e, f) + g
    return f

def gaussian_and_bkg_2(x, a, b, c, d, e):
    "fit function for 356kev peak - cdf fixed to same params as gaussian"
    f = gaussian(x, a, b, c) - d*gaussian_cdf(x, b, c) + e
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

def chi_sq_calc(xdata, ydata, yerr, fit_func, popt):
    "calculate chi sq and p-val of a fit given the data points and fit parameters, e.g. fittype ='linear'"
   
    y_obs = ydata
    y_exp = []

    for index, y_i in enumerate(y_obs):
        x_obs = xdata[index]
        y_exp_i = fit_func(x_obs, *popt)
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

def gauss_count(a,b,c):
    "count/integrate gaussian from -inf to plus inf"
    integral = a*c*np.sqrt(2*np.pi)
    return integral

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
    bguess = xdata[aguess_index] #gauss mean
    cguess =  1 #gauss sigma
    dguess =  min(ydata)/100 #cdf amp
    eguess =  bguess #cdf mean
    fguess = cguess #cdf sigma
    gguess = 0 #offset
    p_guess = [aguess, bguess, cguess, dguess, eguess, fguess, gguess]
    bounds=([0, 0, 0, 0, 0, 0, -np.inf], [np.inf]*7)
    sigma = []
    for index, i in enumerate(yerr):    
        if i != 0:
            sigma.append(yerr[index])
        else:
            sigma.append(1) #just to prevent errors...
    sigma = np.array(sigma)
    popt, pcov = optimize.curve_fit(gaussian_and_bkg, xdata, ydata, p0=p_guess, sigma = sigma, maxfev = 10**7, method ="trf", bounds = bounds)
    
    a,b,c,d,e,f,g = popt[0],popt[1],[2],popt[3],popt[4],popt[5],popt[6]

    fig, ax = plt.subplots()
    #ax.errorbar(xdata, ydata, xerr=0, yerr =yerr, label = "Data", elinewidth = 1, fmt='x', ms = 0.75, mew = 3.0)
    plt.errorbar(xdata, ydata, xerr=0, yerr =yerr, label = "Data", elinewidth = 1, fmt='x', ms = 0.75, mew = 3.0)
    xfit = np.linspace(min(xdata), max(xdata), 1000)
    plt.plot(xfit, gaussian_and_bkg(xfit,*popt), "g", label = "gauss(x,a,b,c) - d*gauss_cdf(x,e,f) + g") 
    plt.plot(xfit, -1*d*gaussian_cdf(xfit,e,f) + g, "r--", label ="-d*gauss_cdf(x,e,f) + g")
    plt.xlabel(key)
    plt.ylabel("Counts")
    plt.legend(loc="upper right", fontsize=8)

    chi_sq, p_value, residuals, dof = chi_sq_calc(xdata, ydata, yerr, gaussian_and_bkg, popt)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    info_str = '\n'.join((r'$a=%.3g \pm %.3g$' % (popt[0], np.sqrt(pcov[0][0])), r'$b=%.3g \pm %.3g$' % (popt[1], np.sqrt(pcov[1][1])), r'$c=%.3g \pm %.3g$' % (popt[2], np.sqrt(pcov[2][2])), r'$d=%.3g \pm %.3g$' % (popt[3], np.sqrt(pcov[3][3])), r'$e=%.3g \pm %.3g$' % (popt[4], np.sqrt(pcov[4][4])), r'$f=%.3g \pm %.3g$' % (popt[5], np.sqrt(pcov[5][5])),r'$g=%.3g \pm %.3g$' % (popt[6], np.sqrt(pcov[6][6])), r'$\chi^2/dof=%.2f/%.0f$'%(chi_sq, dof)))
    plt.text(0.02, 0.98, info_str, transform=ax.transAxes, fontsize=8,verticalalignment='top', bbox=props) #ax.text..ax.tra

    return popt, pcov, xfit


def fit_peak_356_2(key, bins, counts, xmin, xmax):
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
    bguess = xdata[aguess_index] #gauss mean
    cguess =  1 #gauss sigma
    dguess =  min(ydata)/100 #cdf amp
    eguess = 0 #offset
    p_guess = [aguess,bguess,cguess,dguess,eguess]
    bounds=([0, 0, 0, 0, -np.inf], [np.inf]*5)
    #bounds=([0, 0, -np.inf, 0, 0, 0, 0, 0], [np.inf]*4)
    sigma = []
    for index, i in enumerate(yerr):    
        if i != 0:
            sigma.append(yerr[index])
        else:
            sigma.append(1) #just to prevent errors...
    sigma = np.array(sigma)
    popt, pcov = optimize.curve_fit(gaussian_and_bkg_2, xdata, ydata, p0=p_guess, sigma = sigma, maxfev = 10**7, method ="trf", bounds = bounds)
    
    a,b,c,d,e = popt[0],popt[1],[2],popt[3],popt[4]

    fig, ax = plt.subplots()
    #ax.errorbar(xdata, ydata, xerr=0, yerr =yerr, label = "Data", elinewidth = 1, fmt='x', ms = 0.75, mew = 3.0)
    plt.errorbar(xdata, ydata, xerr=0, yerr =yerr, label = "Data", elinewidth = 1, fmt='x', ms = 0.75, mew = 3.0)
    xfit = np.linspace(min(xdata), max(xdata), 1000)
    plt.plot(xfit, gaussian_and_bkg_2(xfit,*popt), "g", label = "gauss(x,a,b,c) - d*gauss_cdf(x,b,c) + e") 
    plt.plot(xfit, -1*d*gaussian_cdf(xfit,b,c) + e, "r--", label ="-d*gauss_cdf(x,b,c) + e")
    plt.xlabel(key)
    plt.ylabel("Counts")
    plt.legend(loc="upper right", fontsize=8)

    chi_sq, p_value, residuals, dof = chi_sq_calc(xdata, ydata, yerr, gaussian_and_bkg_2, popt)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    info_str = '\n'.join((r'$a=%.3g \pm %.3g$' % (popt[0], np.sqrt(pcov[0][0])), r'$b=%.3g \pm %.3g$' % (popt[1], np.sqrt(pcov[1][1])), r'$c=%.3g \pm %.3g$' % (popt[2], np.sqrt(pcov[2][2])), r'$d=%.3g \pm %.3g$' % (popt[3], np.sqrt(pcov[3][3])), r'$e=%.3g \pm %.3g$' % (popt[4], np.sqrt(pcov[4][4])), r'$\chi^2/dof=%.2f/%.0f$'%(chi_sq, dof)))
    plt.text(0.02, 0.98, info_str, transform=ax.transAxes, fontsize=8,verticalalignment='top', bbox=props) #ax.text..ax.tra

    return popt, pcov, xfit


if __name__ =="__main__":
    main()