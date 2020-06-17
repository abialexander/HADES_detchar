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
    #filename = "detector_IC160A_ba_top_01.hdf5"
    filename = "detector_out_am_IC160A_table1_onlygamma_outside_nocollimator_cos_test.hdf5"
    df = pd.read_hdf(filename)
    print(df)


if __name__ =="__main__":
    main()