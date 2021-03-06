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


#this is the OLD SCRIPT - DO NOT USE

def main():

    #read tier 2 runs for Ba data - I02160A
    t2_folder = "/lfs/l1/legend/detector_char/enr/hades/char_data/I02160A/tier2/ba_HS4_top_dlt/pygama/"
    keys, data = read_all_t2(t2_folder)

    no_events = data.size #all events
    print("No. events: ", no_events)

    key = "e_ftp"
    key_data = obtain_key_data(data, keys, key, no_events)

    no_bins = 10000 #7722

    plt.figure()
    counts, bins, bars = plt.hist(key_data, bins=no_bins)
    plt.yscale("log")
    plt.xlabel("e_ftp")
    plt.ylabel("Counts")
    plt.xlim(0, 40000)
    plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/plots/e_ftp.png")

    plt.figure()
    counts, bins, bars = plt.hist(key_data, bins=no_bins)
    plt.yscale("log")
    plt.xlabel("e_ftp")
    plt.ylabel("Counts")
    plt.xlim(0, 7000) 
    plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/plots/e_ftp_zoom.png")


    #___________Calibration__________
    print("Calibrating...")

    print("Calibration peaks:")

    #Fit known calibration peaks from: 
    #https://www.ezag.com/fileadmin/ezag/user-uploads/isotopes/isotopes/Isotrak/isotrak-pdf/Decay_Schema_Data/Ba-133.pdf

    truth_energies = np.array([81.0, 276.40,302.85,356.02,383.85]) #keV
    peak_lims_guess = [[1275, 1315], [4380, 4450], [4805,4870], [5645,5725], [6090,6170]] #rough guess on peak window

    mu_peaks = []
    sigma_peaks = []
    mu_err_peaks = []
    sigma_err_peaks = []

    for index, truth in enumerate(truth_energies):        
        truth_str = str(int(truth))
        print("fitting peak: ", truth_str, " keV")
        xmin, xmax = peak_lims_guess[index][0], peak_lims_guess[index][1]
        print(xmin, " , ", xmax)
        plt.figure()
        mu, sigma, mu_err, sigma_err = fit_peak(key, bins, counts, xmin, xmax)
        mu_peaks.append(mu)
        sigma_peaks.append(sigma)
        mu_err_peaks.append(mu_err)
        sigma_err_peaks.append(sigma_err)
        counts, bins, bars = plt.hist(key_data, bins=no_bins, histtype='step', color='grey')
        plt.xlim(xmin, xmax)
        plt.title(truth_str+" keV peak")
        plt.yscale("log")
        plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/plots/"+truth_str+"keV_peak_uncal.png")


    #Plot calibration curves:
    print("")
    print("Construct calibration curve...")
    calibration_data, calibration_data_err = np.array(mu_peaks), np.array(mu_err_peaks)
    plt.figure()
    m, c, m_err, c_err, chi_sq, p_value, residuals, dof = calibration(calibration_data, calibration_data_err, truth_energies, fittype="linear_fit")
    for x, y in zip(truth_energies, calibration_data):
        truth_str = str(int(x))
        plt.annotate(truth_str, (x,y), textcoords="offset points", xytext=(-5,5), ha='center') # horizontal alignment can be left, right or center
    plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/plots/calibration_curve_linear.png")

    #plot residuals
    plt.figure()
    residuals_err = [1.0]*len(residuals) #error on residuals is 1 by definition
    plt.errorbar(truth_energies, residuals, yerr = residuals_err, fmt ='bo') #, ms = 1.0, mew = 3.0)
    plt.ylabel("r=(data-fit)/error")
    plt.xlabel("Energy (keV)")
    plt.title("Residual Plot for Calibration Graph")
    for x, y in zip(truth_energies, residuals):
        truth_str = str(int(x))
        plt.annotate(truth_str, (x,y), textcoords="offset points", xytext=(0,5), ha='center') # horizontal alignment can be left, right or center 
    plt.axhline(linewidth=2, color='black', dashes = (5,2,1,2))
    plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/plots/calibration_residuals.png")


    #try non-linear calibration
    plt.figure()
    a, b, c_quad, a_err, b_err, c_quad_err, chi_sq, p_value, residuals, dof = calibration(calibration_data, calibration_data_err, truth_energies, fittype="quadratic_fit")
    for x, y in zip(truth_energies, calibration_data):
        truth_str = str(int(x))
        plt.annotate(truth_str, (x,y), textcoords="offset points", xytext=(-5,5), ha='center') 
    plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/plots/calibration_curve_quadratic.png")



    #Save calibration coefficients to json file
    calibration_coef_dict = {
        "m": m,
        "m_err" : m_err,
        "c" : c,
        "c_err" : c_err
    }
    with open("calibration_coef.json", "w") as outfile: 
        json.dump(calibration_coef_dict, outfile)

    #Linearly calibrated data:
    print("")
    print("Linearly calibrate energy...")

    calibrated_energy = (key_data-c)/m

    bin_width = 0.5 #0.5 kev = resolution
    no_bins_ideal = int(max(calibrated_energy/bin_width))
    #print("ideal no bins: ", no_bins_ideal) #=7722

    plt.figure()
    counts, bins_cal, bars = plt.hist(calibrated_energy, bins=no_bins)
    plt.xlabel("Energy (KeV)")
    plt.ylabel("Frequency")
    plt.xlim(0, 2500)
    plt.yscale("log")
    plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/plots/calibrated_energy.png") 

    #plot zoomed in
    plt.figure()
    counts, bins_cal, bars = plt.hist(calibrated_energy, bins=100000)
    plt.xlabel("Energy (KeV)")
    plt.ylabel("Frequency")
    plt.yscale("log")
    plt.xlim(0,450)
    plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/plots/calibrated_energy_zoom.png") 
    
    #__________CALIBRATED___________




    # #Calibrate - use linear calibration coefficients from previous work
    # m, m_err = 15.786, 0.007
    # c, c_err = -387.802, 11.591



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

    #print(df_total.describe())
    #print(df_total.head())

    # plt.figure()
    # df_total.hist('e_ftp',  bins = 10000) #testing
    # plt.ylabel("Frequency")
    # plt.yscale("log")
    # plt.savefig("/lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/plots/e_ftp_pandas.png")

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

# def gaussian_cdf_bkg(x,a,b,c,d):

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
    f = a*np.sqrt(x) +c
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


def calibration(data, data_err, truth, fittype):
    "calibrate energy data for given peaks, fittype = linear_fit or quadratic_fit"
    xdata = truth
    ydata = data
    yerr = data_err

    if fittype == "linear_fit":
        popt, pcov = optimize.curve_fit(linear_fit, xdata, ydata, sigma = yerr, maxfev = 1000000, method ="trf")
        m, c = popt[0], popt[1]
        m_err, c_err = np.sqrt(pcov[0][0]), np.sqrt(pcov[1][1])
        fig, ax = plt.subplots()
        xfit = np.linspace(min(xdata), max(xdata), 1000)
        yfit = linear_fit(xfit,*popt)
        plt.plot(xfit, yfit, "g", label = "linear fit")
        ax.errorbar(truth, data, xerr=0, yerr =yerr, elinewidth = 1, fmt='x', ms = 1.5, mew = 4.0, label = "callibration peaks")
        plt.xlabel("truth (keV)")
        plt.ylabel("calibration peaks")
        plt.legend(loc = 4)

        chi_sq, p_value, residuals, dof = chi_sq_calc(xdata, ydata, yerr, "linear", popt)

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        info_str = '\n'.join((r'$m = %.3f \pm %.3f$' % (m, m_err, ), r'$c = %.3f \pm %.3f$' % (c, c_err,), r'$\chi^2/dof=%.0f/%.0f$'%(chi_sq, dof)))# r'$\chi^2 = %.3f$'%(chi_sq), r'$p = %.3g$'%(p_value)))
        ax.text(0.05, 0.95, info_str, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)

        return m, c, m_err, c_err , chi_sq, p_value, residuals, dof

    if fittype == "quadratic_fit":
        popt, pcov = optimize.curve_fit(quadratic_fit, xdata, ydata, sigma = yerr, maxfev = 1000000, method ="trf")
        a,b,c = popt[0], popt[1], popt[2]
        a_err, b_err, c_err = np.sqrt(pcov[0][0]), np.sqrt(pcov[1][1]), np.sqrt(pcov[2][2])
        fig, ax = plt.subplots()
        xfit = np.linspace(min(xdata), max(xdata), 1000)
        yfit = quadratic_fit(xfit,*popt)
        plt.plot(xfit, yfit, "g", label = "quadratic fit")
        ax.errorbar(truth, data, xerr=0, yerr =yerr, elinewidth = 1, fmt='x', ms = 1.5, mew = 4.0, label = "callibration peaks")
        plt.xlabel("truth (keV)")
        plt.ylabel("calibration peaks")
        plt.legend(loc = 4)

        chi_sq, p_value, residuals, dof = chi_sq_calc(xdata, ydata, yerr, "quadratic", popt)

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        info_str = '\n'.join((r'$a = %.3g \pm %.3g$' % (a, a_err, ), r'$b = %.3f \pm %.3f$' % (b, b_err,), r'$c = %.3f \pm %.3f$' % (c, c_err,),  r'$\chi^2/dof=%.0f/%.0f$'%(chi_sq, dof))) #, r'$p = %.3g$'%(p_value)))
        ax.text(0.05, 0.95, info_str, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)

        return a, b, c, a_err, b_err, c_err, chi_sq, p_value, residuals, dof


if __name__ =="__main__":
    main()