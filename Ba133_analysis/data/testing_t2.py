import pandas as pd
import numpy as np
import h5py
import os
import matplotlib.pyplot as plt
import argparse
from scipy import optimize
from scipy import stats
import glob

import seaborn as sb
from sklearn.model_selection import train_test_split 
from sklearn.linear_model import LinearRegression
from sklearn import metrics

#import ROOT as root
"Analyses 2 files concatonated together as a test"

def main():

    #Get data
    key = "e_ftp"

    #test on 2 files
    #keys, data = read_2_t2()
    #keys_t1, data_t1 = read_2_t1()

    #or whole data set
    keys, data = read_all_t2()
    #keys_t1, data_t1 = read_all_t1()

    no_events = data.size #all events
    key_data = obtain_key_data(data, keys, key, no_events)
    plt.figure()
    counts, bins, bars = plot_histo_matplot(key_data, key)
    plt.yscale("log")
    #plt.savefig("/unix/legend/abi/I02160A/tier2/th_HS2_top_psa/pygama/plots/"+key+".png") 


    #_________END POINT_____ - tl-208 peak (2614)
    print("")
    print("Tl-208 Endpoint: ")
    xmin_ep, xmax_ep = 40800., 41000. #by inspection
    mu_ep, sigma_ep, mu_err_ep, sigma_err_ep = fit_peak(key, bins, counts, xmin_ep, xmax_ep)
    plt.title("Tl-208 End point (2614 KeV)")
    #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/EP_test.png")
    plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/EP.png")

    # peak_ievt_ep = identify_peak_ievt(mu_ep, sigma_ep, key_data, data, key, keys) #identify peak wfs

    # #ep_df = pd.DataFrame(peak_ievt_ep)
    # #ep_df.to_csv("ep_wf_list_test.csv", sep=",", index=False)

    # plt.figure()
    # wf_list, wf_max_list = plot_peak_wfs(peak_ievt_ep, data_t1, keys_t1, "EP")
    # plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/av_wf_EP_test.png")
    # #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/av_wf_EP.png")

    # #wf max - testing
    # plt.figure()
    # plt.hist(wf_max_list, bins = 1000)
    # plt.xlabel("wf max")
    # plt.ylabel("frequency")
    # plt.title("Wf maxima - End Point")
    # plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/av_wf_max_EP_test.png")
    # #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/av_wf_max_EP_test1run.png")
    # #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/av_wf_max_EP.png")

    #________SINGLE ESCAPE PEAK____________ (2103)
    print("")
    print("Tl-208 SEP: ")
    xmin_sep, xmax_sep = 32725., 32900.
    mu_sep, sigma_sep, mu_err_sep, sigma_err_sep = fit_peak(key, bins, counts, xmin_sep, xmax_sep)
    plt.title("Tl-208 Single escape peak (2103 KeV)")
    #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/SEP_test.png")
    plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/SEP.png")

    # peak_ievt_sep = identify_peak_ievt(mu_sep, sigma_sep, key_data, data, key, keys) #identify peak wfs
    # plt.figure()
    # wf_list, wf_max_list = plot_peak_wfs(peak_ievt_sep, data_t1, keys_t1, "SEP")
    # plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/av_wf_SEP_test.png")
    # #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/av_wf_SEP.png")

    # plt.figure()
    # plt.hist(wf_max_list, bins = 1000)
    # plt.xlabel("wf max")
    # plt.ylabel("frequency")
    # plt.title("Wf maxima - SEP")
    # plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/av_wf_max_SEP_test.png")
    # #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/av_wf_max_SEP_test1run.png")
    # #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/av_wf_max_SEP.png")

    #___________BI-212a____________ (1620)
    print("")
    print("Bi-212 a: ")
    xmin_Bi212a, xmax_Bi212a = 25125., 25250.
    mu_Bi212a, sigma_Bi212a, mu_err_Bi212a, sigma_err_Bi212a = fit_peak(key, bins, counts, xmin_Bi212a, xmax_Bi212a)
    plt.title("Bi-212 a (1620 KeV)")
    #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/Bi212a_test.png")
    plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/Bi212a.png")


    #_________Double Escape Peak________ (1592)
    print("")
    print("Tl-208 DEP: ")
    xmin_dep, xmax_dep = 24675., 24800.
    mu_dep, sigma_dep, mu_err_dep, sigma_err_dep = fit_peak(key, bins, counts, xmin_dep, xmax_dep)
    plt.title("Tl-208 Double Escape Peak (1592 KeV)")
    #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/DEP_test.png")
    plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/DEP.png")
    
    # peak_ievt_dep = identify_peak_ievt(mu_dep, sigma_dep, key_data, data, key, keys) #identify peak wfs
    # plt.figure()
    # wf_list, wf_max_list = plot_peak_wfs(peak_ievt_dep, data_t1, keys_t1, "DEP")
    # plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/av_wf_DEP_test.png")
    # #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/av_wf_DEP.png")

    # plt.figure()
    # plt.hist(wf_max_list, bins = 1000)
    # plt.xlabel("wf max")
    # plt.ylabel("frequency")
    # plt.title("Wf maxima - DEP")
    # plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/av_wf_max_DEP_test.png")
    # #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/av_wf_max_DEP_test1run.png")
    # #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/av_wf_max_DEP.png")

    #__________BI 212b_ ______ (1512)
    print("")
    print("Bi-212: ")
    xmin_Bi212b, xmax_Bi212b = 23425., 23525.
    mu_Bi212b, sigma_Bi212b, mu_err_Bi212b, sigma_err_Bi212b = fit_peak(key, bins, counts, xmin_Bi212b, xmax_Bi212b)
    plt.title("Bi-212 b (1512 KeV)")
    #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/Bi212b_test.png")
    plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/Bi212b.png")


    #__________Combined WFs for Tl-208______
    # plt.figure() #combined wfs
    # plot_peak_wfs(peak_ievt_ep, data_t1, keys_t1, "EP")
    # plot_peak_wfs(peak_ievt_sep, data_t1, keys_t1, "SEP")
    # plot_peak_wfs(peak_ievt_dep, data_t1, keys_t1, "DEP")
    # plt.legend()
    # plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/av_wf_all_test.png")
    # #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/av_wf_all_test1run.png")
    # #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/av_wf_all.png")

    # plt.figure() #combined wfs -  zoomed in
    # plot_peak_wfs(peak_ievt_ep, data_t1, keys_t1, "EP")
    # plot_peak_wfs(peak_ievt_sep, data_t1, keys_t1, "SEP")
    # plot_peak_wfs(peak_ievt_dep, data_t1, keys_t1, "DEP")
    # plt.legend()
    # #plt.xlim(1750,1950) 
    # plt.xlim(25,35) #- micros
    # plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/av_wf_all_zoom_test.png")
    # #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/av_wf_all_zoom.png")


    plt.close('all')

    # #________-CALIBRATION__________
    print("")
    print("Calibration: ")

    #linear calibration
    m_e = 510.99 #electron mass, kev
    calibration_data = np.array([mu_ep, mu_sep, mu_Bi212a, mu_dep, mu_Bi212b])
    calibration_data_err = np.array([mu_err_ep, mu_err_sep, mu_err_Bi212a, mu_err_dep, mu_err_Bi212b])
    truth = np.array([2614.511, 2614.511-m_e, 1620, 2614.511-2*m_e, 1512]) #actual energy values for these peaks, https://www.nndc.bnl.gov/nudat2/decaysearchdirect.jsp?nuc=208TL&unc=nds
    
    #using scipy.curvefit
    plt.figure()
    m, c, m_err, c_err, chi_sq, p_value, residuals = calibration(calibration_data, calibration_data_err, truth, fittype="linear_fit")
    labels = ["Tl-208", "SEP", "Bi-212", "DEP", "Bi-212"]
    for x, y in zip(truth, calibration_data):
        index = np.where(truth == x)[0][0]
        label = labels[index]
        plt.annotate(label, (x,y), textcoords="offset points", xytext=(-5,5), ha='center') # horizontal alignment can be left, right or center
    #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/calibration_test.png")
    plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/calibration.png")

    #plot residuals
    plt.figure()
    residuals_err = [1.0]*len(residuals) #error on residuals is 1 by definition
    plt.errorbar(truth, residuals, yerr = residuals_err, fmt ='bo')
    plt.ylabel("r=(data-fit)/error")
    plt.xlabel("Energy (keV)")
    plt.title("Residual Plot for Calibration Graph")
    for x, y in zip(truth, residuals):
        index = np.where(truth == x)[0][0]
        label = labels[index]
        plt.annotate(label, (x,y), textcoords="offset points", xytext=(0,5), ha='center') # horizontal alignment can be left, right or center 
    plt.axhline(linewidth=2, color='black', dashes = (5,2,1,2))
    #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/calibration_residuals_test.png")
    plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/calibration_residuals.png")

    #try non-linear calibration
    plt.figure()
    a, b, c_quad, a_err, b_err, c_quad_err, chi_sq, p_value, residuals = calibration(calibration_data, calibration_data_err, truth, fittype="quadratic_fit")
    labels = ["Tl-208", "SEP", "Bi-212", "DEP", "Bi-212"]
    for x, y in zip(truth, calibration_data):
        index = np.where(truth == x)[0][0]
        label = labels[index]
        plt.annotate(label, (x,y), textcoords="offset points", xytext=(-5,5), ha='center') # horizontal alignment can be left, right or center
    #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/calibration_quadratic_test.png")
    plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/calibration_quadratic.png")


    #using seaborn linear regress -testing
    # plt.figure()
    # calibration(calibration_data, calibration_data_err, truth, linreg = True)
    # labels = ["Tl-208", "SEP", "Bi-212", "DEP", "Bi-212"]
    # for x, y in zip(truth, calibration_data):
    #     index = np.where(truth == x)[0][0]
    #     label = labels[index]
    #     plt.annotate(label, (x,y), textcoords="offset points", xytext=(-5,5), ha='center') # horizontal alignment can be left, right or center
    # #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/calibration_test_linreg.png")


    #___________Calibrated Graphs: LINEAR____________
    print("")
    print("Calibrated graphs: ")

    #calibrated energy spectrum 
    calibrated_energy = (key_data-c)/m
    plt.figure()
    counts, bins_cal, bars = plot_histo_matplot(calibrated_energy, "Energy (KeV)")
    plt.yscale("log")
    #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/calibrated_energy_spectrum_test.png") 
    plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/calibrated_energy_spectrum.png") 


    #calibrated peaks
    key = "Energy (KeV)"

    #tl-end point peak
    print("")
    print("Tl end point: ")
    xmin_ep_cal, xmax_ep_cal = (xmin_ep-c)/m, (xmax_ep-c)/m
    mu_ep_cal, sigma_ep_cal, mu_err_ep_cal, sigma_err_ep_cal = fit_peak(key, bins_cal, counts, xmin_ep_cal, xmax_ep_cal)
    plt.title("Tl-208 End point (2614 KeV) - Calibrated")
    #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/EP_cal_test.png")
    plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/EP_cal.png")

    #tl sep
    print("")
    print("Tl SEP: ")
    xmin_sep_cal, xmax_sep_cal = (xmin_sep-c)/m, (xmax_sep-c)/m
    mu_sep_cal, sigma_sep_cal, mu_err_sep_cal, sigma_err_sep_cal = fit_peak(key, bins_cal, counts, xmin_sep_cal, xmax_sep_cal)
    plt.title("Tl-208 Single Escape Peak (2103 KeV) - Calibrated")
    #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/SEP_cal_test.png")
    plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/SEP_cal.png")

    #tl- dep
    print("")
    print("Tl dep: ")
    xmin_dep_cal, xmax_dep_cal = (xmin_dep-c)/m, (xmax_dep-c)/m
    mu_dep_cal, sigma_dep_cal, mu_err_dep_cal, sigma_err_dep_cal = fit_peak(key, bins_cal, counts, xmin_dep_cal, xmax_dep_cal)
    plt.title("Tl-208 Double Escape Peak (1592 KeV) - Calibrated")
    #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/DEP_cal_test.png")
    plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/DEP_cal.png")

    #Bi212 a
    print("")
    print("Bi-212 a: ")
    xmin_Bi212a_cal, xmax_Bi212a_cal = (xmin_Bi212a-c)/m, (xmax_Bi212a-c)/m
    mu_Bi212a_cal, sigma_Bi212a_cal, mu_err_Bi212a_cal, sigma_err_Bi212a_cal = fit_peak(key, bins_cal, counts, xmin_Bi212a_cal, xmax_Bi212a_cal)
    plt.title("Bi-212 a (1620 KeV) - Calibrated")
    #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/Bi212a_cal_test.png")
    plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/Bi212a_cal.png")

    #Bi212 b
    print("")
    print("Bi-212 b: ")
    xmin_Bi212b_cal, xmax_Bi212b_cal = (xmin_Bi212b-c)/m, (xmax_Bi212b-c)/m
    mu_Bi212b_cal, sigma_Bi212b_cal, mu_err_Bi212b_cal, sigma_err_Bi212b_cal = fit_peak(key, bins_cal, counts, xmin_Bi212b_cal, xmax_Bi212b_cal)
    plt.title("Bi-212 b (1512 KeV) - Calibrated")
    #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/Bi212b_cal_test.png")
    plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/Bi212b_cal.png")


    #________RESOLUTION___________
    #plot FWHM of calibration peaks against energy
    print("")
    print("Resolution: ")

    #This method does not work...
    # energies = np.array(truth)
    # sigmas_uncal = [sigma_ep, sigma_sep, sigma_Bi212a, sigma_dep, sigma_Bi212b]
    # sigmas = (np.array(sigmas_uncal)-c)/m
    # sigmas_err_uncal = [sigma_err_ep, sigma_err_sep, sigma_err_Bi212a, sigma_err_dep, sigma_err_Bi212b]
    # sigmas_err = (np.array(sigmas_err_uncal)-c)/m
    # FWHM = 2*np.sqrt(2*np.log(2))*(np.array(sigmas)) #~2.35*sigma
    # FWHM_err = 2*np.sqrt(2*np.log(2))*(np.array(sigmas_err))

    # plt.figure()
    # plt.errorbar(energies, FWHM, xerr=0, yerr =FWHM_err, label = "Data", fmt = 'o')
    # plt.xlabel("Energy (KeV)")
    # plt.ylabel("FWHM (KeV)")
    # plt.legend()
    # for x, y in zip(energies, FWHM):
    #     index = np.where(energies == x)[0][0]
    #     label = labels[index]
    #     plt.annotate(label, (x,y), textcoords="offset points", xytext=(-10,5), ha='center') # horizontal alignment can be left, right or center
    # plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/resolution_test.png")
    # #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/resolution.png")


    energies = [mu_ep_cal, mu_sep_cal, mu_Bi212a_cal, mu_dep_cal, mu_Bi212b_cal]
    energies_err = [mu_err_ep_cal, mu_err_sep_cal, mu_err_Bi212a_cal, mu_err_dep_cal, mu_err_Bi212b_cal]
    sigmas = np.array([sigma_ep_cal, sigma_sep_cal, sigma_Bi212a_cal, sigma_dep_cal, sigma_Bi212b_cal])
    sigmas_err = np.array([sigma_err_ep_cal, sigma_err_sep_cal, sigma_err_Bi212a_cal, sigma_err_dep_cal, sigma_err_Bi212b_cal])
    FWHM = 2*np.sqrt(2*np.log(2))*(np.array(sigmas)) #~2.35*sigma
    FWHM_err = 2*np.sqrt(2*np.log(2))*(np.array(sigmas_err))

    #sqrt fit - excluding sep
    xdata, ydata, yerr = energies, FWHM, FWHM_err
    np.delete(xdata,1)
    np.delete(ydata,1) #remove SEP from fit
    np.delete(yerr, 1)
    Aguess = max(ydata) - min(ydata)
    offset_guess = 0
    p_guess = [Aguess, offset_guess]
    bounds=(0, [np.inf, np.inf])
    sigma = []
    for index, i in enumerate(yerr):    
        if i != 0:
            sigma.append(yerr[index])
        else:
            sigma.append(1) #just to prevent errors...
    sigma = np.array(sigma)
    popt, pcov = optimize.curve_fit(sqrt_curve, xdata, ydata, p0=p_guess, sigma = sigma, maxfev = 1000000, method ="trf", bounds = bounds)
    A, offset = popt[0], popt[1] #must be positive
    A_err, offset_err = np.sqrt(pcov[0][0]), np.sqrt(pcov[1][1])

    fig, ax = plt.subplots()
    ax.errorbar(energies, FWHM, xerr=energies_err, yerr = FWHM_err, label = "Data", fmt='o') #, elinewidth = 1, fmt='x', ms = 1.5, mew = 4.0)
    xfit = np.linspace(min(xdata), max(xdata), 1000)
    plt.plot(xfit, sqrt_curve(xfit,*popt), "g", label = "$A*\sqrt{x} + c$ fit")
    plt.xlabel("Energy (keV)")
    plt.ylabel("FWHM (keV)")
    plt.legend()
    chi_sq, p_value, residuals = chi_sq_calc(xdata, ydata, yerr, "sqrt", popt)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    info_str = '\n'.join((r'$A=%.3f \pm %.3f$' % (A, A_err, ), r'$c=%.3f \pm %.3f$' % (offset, offset_err,), r'$\chi^2=%.3f$'%chi_sq, r'$p=%.3g$'%p_value))
    ax.text(0.05, 0.95, info_str, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)
    for x, y in zip(energies, FWHM):
        index = np.where(energies == x)[0][0]
        label = labels[index]
        plt.annotate(label, (x,y), textcoords="offset points", xytext=(-10,5), ha='center') # horizontal alignment can be left, right or center
   
    #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/resolution_cal_test.png")
    plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/resolution_cal.png")



    #plt.show()

def read_t2(filepath):
    "get data from t2 hdf5 file"

    df = pd.read_hdf(filepath, "data") 
    df.reset_index(inplace=True)
    keys = df.keys()
    print("Available keys are: ")
    print(keys)
    data = df.to_numpy()
    print("data shape: ", data.shape)

    return keys, data

def read_2_t2():
    "for testing only"

    directorypath = "/unix/legend/abi/I02160A/tier2/th_HS2_top_psa/pygama/"
    file1 = "t2_char_data-I02160A-th_HS2_top_psa-run0001-191021T162944.h5"
    file2 = "t2_char_data-I02160A-th_HS2_top_psa-run0001-191021T163144.h5"

    df1 = pd.read_hdf(directorypath+file1, "data")
    df1.reset_index(inplace=True)
    df2 = pd.read_hdf(directorypath+file2, "data")
    df2.reset_index(inplace=True)
    df_total = pd.concat([df1, df2], axis=0, ignore_index=True)
    df_total.reset_index(inplace=True)

    #print(df_total.describe())
    #print(df_total.head())
    # plt.figure()
    # df1.hist('e_ftp',  bins = 10000) #testing
    # plt.ylabel("Frequency")
    #plt.savefig("/home/aalexander/LEGEND/HADES_analysis/th228_top_plots/e_ftp_pandas_test.png")

    #df_total = df1 #testing with 1 file
    keys = df_total.keys()
    print("Available keys are: ")
    print(keys)
    data = df_total.to_numpy()
    print("data shape: ", data.shape)

    return keys, data

def read_all_t2():
    "get data from all t2 files from same run within a directory"

    directorypath_t2 = "/unix/legend/abi/I02160A/tier2/th_HS2_top_psa/pygama/"

    run1_files = glob.glob(directorypath_t2 + "/*run0001*.h5")

    file_list = []
    for filename in run1_files:
        df = pd.read_hdf(filename, "data")
        file_list.append(df)
    
    df_total = pd.concat(file_list, axis=0, ignore_index=True)
    keys = df_total.keys()
    data = df_total.to_numpy()

    return keys, data

def read_2_t1():
    "get data from 2 t1 hdf5 file"

    directorypath = "/unix/legend/abi/I02160A/tier1/th_HS2_top_psa/pygama/"
    file1 = "t1_char_data-I02160A-th_HS2_top_psa-run0001-191021T162944.h5"
    file2 = "t1_char_data-I02160A-th_HS2_top_psa-run0001-191021T163144.h5"

    df1 = pd.read_hdf(directorypath+file1)
    df1.reset_index(inplace=True)
    df2 = pd.read_hdf(directorypath+file2)
    df2.reset_index(inplace=True)
    df_total = pd.concat([df1, df2], axis=0, ignore_index=False)
    df_total.reset_index(inplace=True)

    #df_total = df1 #test with just 1 file
    data = df_total.to_numpy()
    keys = df_total.keys()

    return keys, data

def read_all_t1():
    "get data from all t1 files from same run within a directory"

    directorypath_t1 = "/unix/legend/abi/I02160A/tier1/th_HS2_top_psa/pygama/"
    
    run1_files = glob.glob(directorypath_t1 + "/*run0001*.h5")

    file_list = []
    for filename in run1_files:
        df = pd.read_hdf(filename)
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


def plot_histo_root(key_data, key):
    "plot histo of from a list of a quantity using root - needs pyroot"

    #list = ["5","9"]
    hist = root.TH1D("hist name","hist", 100, 0, 100)
    for i in range(0,hist.GetNbinsX()):
        value = key_data[i]
        hist.SetBinContent(i+1,value)
    
    canvas = root.TCanvas( 'c1', 'Histogram', 200, 10, 700, 900 )
    canvas.cd ()
    pad = root.TPad( 'pad1', 'The pad with the function',  0.03, 0.62, 0.50, 0.92, 21 )
    pad.Draw()
    hist.Draw('h')

    

def plot_histo_matplot(key_data, key):
    "plot histo from list of a quantity using matplotlib"

    #plt.figure()
    counts, bins, bars = plt.hist(key_data, bins=10000)
    plt.xlabel(key)
    plt.ylabel("Frequency")

    return counts, bins, bars


def gaussian(x,a,b,c,d):
    "gaussian function with offset d"
    f = a*np.exp(-((x-b)**2.0/(2.0*c**2.0))) +d
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
    deg_freedom = N-1
    chi_sq_red = chi_sq/deg_freedom

    p_value = 1-stats.chi2.cdf(chi_sq, deg_freedom)
    

    print("chisq with errors:")
    print("chisq: ", chi_sq)
    print("dof: ", deg_freedom)
    print("p: ", p_value)

    # chi_sq, p_value = stats.chisquare(y_obs, y_exp)#this is the old chi_sq without errors

    return chi_sq, p_value, residuals

def fit_peak(key, bins, counts, xmin, xmax): #p_guess):
    "fit a gaussian to a peak and return mean and FWHM range"

    #counts, bins, bars = plt.hist(key_data, bins=10000)
    no_bins = bins.size #=initially 10000

    xdata = []
    ydata = []
    for bin in bins:
        #bin = bin - 0.5*(max(bins)-(min(bins)))/no_bins #this leads to incorrect indexing so leave out for now
        if bin < xmax and bin > xmin:
            xdata.append(bin)
            bin_index = np.where(bins == bin)[0][0]
            ydata.append(counts[bin_index])

    xdata = np.array(xdata)
    ydata = np.array(ydata)     
    #xdata = xdata.astype(float)
    #print(xdata.dtype)

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
    ax.errorbar(xdata, ydata, xerr=0, yerr =yerr, label = "Data", elinewidth = 1, fmt='x', ms = 1.5, mew = 4.0)
    xfit = np.linspace(min(xdata), max(xdata), 1000)
    plt.plot(xfit, gaussian(xfit,*popt), "g", label = "Gaussian fit")
    plt.xlabel(key)
    plt.ylabel("Frequency")
    plt.legend()

    chi_sq, p_value, residuals = chi_sq_calc(xdata, ydata, yerr, "gaussian", popt)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    info_str = '\n'.join((r'$\mu=%.3f \pm %.3f$' % (mu, mu_err, ), r'$\sigma=%.3f \pm %.3f$' % (np.abs(sigma), sigma_err,), r'$\chi^2=%.3f$'%chi_sq, r'$p=%.3g$'%p_value))
    ax.text(0.05, 0.95, info_str, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)

    FWHM = 2*np.sqrt(2*np.log(2))*sigma #gaussian FWHM relationship
    peak_range = [mu-FWHM, mu+FWHM]

    return mu, sigma, mu_err, sigma_err

def identify_peak_ievt(mu, sigma, key_data, data, key, keys):
    "find event number of events within the FWHM of a gaussian peak"

    #find energies of events within peak range: mu +/- FWHM
    FWHM = 2*np.sqrt(2*np.log(2))*sigma #gaussian FWHM relationship
    peak_range = [mu-FWHM, mu+FWHM]

    peak_energies = []
    for i in key_data:
        #if i > mu-FWHM and i < mu+FWHM:
        if i > mu-sigma and i<mu+sigma:
        #if i > mu-0.1*sigma and i<mu+0.1*sigma: #tighter restriction
            peak_energies.append(i)
    
    #find event number of events within this range
    peak_ievt = []
    key_index = (np.where(keys == key))[0][0]
    ievt_index = (np.where(keys == "ievt"))[0][0]
    for i in peak_energies:
        energy_index = (np.where(key_data == i))[0][0]
        event_data = data[energy_index,:]
        ievt = event_data[ievt_index]
        peak_ievt.append(ievt)

    ievt_list = data[:, ievt_index]
    return peak_ievt

def plot_peak_wfs(peak_ievt, data_t1, keys_t1, label):

    index0 = (np.where(keys_t1 == 0))[0][0]
    index3747 = (np.where(keys_t1 == 3747))[0][0]

    time = keys_t1[index0:index3747]

    ievt_key_index = np.where(keys_t1 == "ievt")[0][0]
    ievt_list = data_t1[:,ievt_key_index]
    wf_list = []
    wf_max_list = []
    wf_nos = []
    for ievt in peak_ievt:
        #wf_no = (np.where(ievt_list == ievt))[0][0]
        wf_no = int(ievt_list[ievt])
        wf_nos.append(wf_no)
        wf = data_t1[wf_no][index0:index3747]
        wf_max = max(wf)
        wf_list.append(wf)
        wf_max_list.append(wf_max)
    
    wf_list = np.array(wf_list)
    wf_max_list = np.array(wf_max_list)
    #average peak wfs:
    av_wf = np.mean(wf_list, axis=0)
    time = (np.array(time)*16)/1000 #sample every 16 ns
    plt.plot(time, av_wf, label = label)
    plt.xlabel("time ($\mu$s)")
    plt.ylabel("Average WF")

    return wf_list, wf_max_list

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

        chi_sq, p_value, residuals = chi_sq_calc(xdata, ydata, yerr, "linear", popt)

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        info_str = '\n'.join((r'$m = %.3f \pm %.3f$' % (m, m_err, ), r'$c = %.3f \pm %.3f$' % (c, c_err,), r'$\chi^2 = %.3f$'%(chi_sq), r'$p = %.3g$'%(p_value)))
        ax.text(0.05, 0.95, info_str, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)

        return m, c, m_err, c_err , chi_sq, p_value, residuals

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

        chi_sq, p_value, residuals = chi_sq_calc(xdata, ydata, yerr, "quadratic", popt)

        print("a: ", a)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        info_str = '\n'.join((r'$a = %.3g \pm %.3g$' % (a, a_err, ), r'$b = %.3f \pm %.3f$' % (b, b_err,), r'$c = %.3f \pm %.3f$' % (c, c_err,), r'$\chi^2 = %.3f$'%(chi_sq), r'$p = %.3g$'%(p_value)))
        ax.text(0.05, 0.95, info_str, transform=ax.transAxes, fontsize=10,verticalalignment='top', bbox=props)

        return a, b, c, a_err, b_err, c_err, chi_sq, p_value, residuals


    #linear regression method - test
    #sb.set(color_codes=True)
    #ax = sb.regplot(x=xdata, y=ydata)
    #plt.xlabel("truth (keV)")
    #plt.ylabel("calibration peaks")
    #plt.legend(loc = 4)



if __name__ =="__main__":
    main()