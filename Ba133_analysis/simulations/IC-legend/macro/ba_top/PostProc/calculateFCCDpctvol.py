import sys, h5py
import pandas as pd
import numpy as np
from datetime import datetime


def main():  

    #print date and time for log:
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S") # dd/mm/YY H:M:S
    print("")
    print("date and time =", dt_string)	
    print("")


    FCCD = 0.57 #mm

    #define geometry
    position_crystal_from_top = 7.0 #all in mm, tajen from gdml files
    crystal_height = 65.4
    bottom_height = 45.3
    cavity_height = 33.7
    groove_height = 2.0
    top_rad = 75.5/2 #top_width/2
    bottom_rad = 79.8/2
    cavity_rad = 9.3/2
    groove_inner_width = 22.6 #ignore grooves for now
    groove_outer_width = 31.0


    
    #CALCULATE DET VOLUME

    FCCD = 0. 
    #region 1 geometry
    r = bottom_rad - FCCD
    l = crystal_height - bottom_height - FCCD
    H = position_crystal_from_top + FCCD
    A = top_rad - FCCD
    B = bottom_rad - FCCD
    z0 = H - l*A/(B-A)
    h = l + H - z0

    r_cylinder = cavity_rad + FCCD

    #region 2 geometry
    l2 = cavity_height - (l+H)

    #region 3 geometry
    l3 = bottom_height - l2 - FCCD


    detVol_1 = (1/3)*np.pi*bottom_rad**2*h - (1/3)*np.pi*top_rad**2*(h-l)
    detVol_2 = np.pi*bottom_rad**2*bottom_height
    detVol = detVol_1 + detVol_2 - np.pi*cavity_rad**2*cavity_height
    print(detVol, "mm^3")
    print(detVol/10**9, "m^3")


    #calculate FCCD volume
    FCCD = 0.568 
    #FCCD_vol_1 = np.pi*()*l
    r1, r2 = bottom_rad, bottom_rad-FCCD
    h1, h2, hhat = h, h-FCCD, l
    FCCD_vol_1 = (1/3)*np.pi*(r1**2*h1 - r2**2*h2 - (h1-hhat)*((h1-hhat)*(r1/h1))**2 +(h2-hhat)*((h2-hhat)*(r2/h2))**2)
    FCCD_vol_2 = np.pi*(bottom_rad**2-(bottom_rad-FCCD)**2)*(bottom_height-FCCD) 
    FCCD_vol = FCCD_vol_1 + FCCD_vol_2 + np.pi*((cavity_rad+0.5*FCCD)**2-cavity_rad**2) + np.pi*(cavity_rad+FCCD)**2*(0.5*FCCD)

    print(FCCD_vol, "mm^3")
    print(FCCD_vol/10**9, "m^3")

    FCCDpctvol =FCCD_vol/detVol
    print(FCCDpctvol*100, " %")



    #calculate error pct volume
    FCCD_err = 0.046
    r1, r2 = bottom_rad, bottom_rad-FCCD_err
    h1, h2, hhat = h, h-FCCD_err, l
    FCCD_vol_1_err = (1/3)*np.pi*(r1**2*h1 - r2**2*h2 - (h1-hhat)*((h1-hhat)*(r1/h1))**2 +(h2-hhat)*((h2-hhat)*(r2/h2))**2)
    FCCD_vol_2_err = np.pi*(bottom_rad**2-(bottom_rad-FCCD_err)**2)*(bottom_height-FCCD_err) 
    FCCD_vol_err = FCCD_vol_1_err + FCCD_vol_2_err + np.pi*((cavity_rad+0.5*FCCD_err)**2-cavity_rad**2) + np.pi*(cavity_rad+FCCD_err)**2*(0.5*FCCD_err)

    print(FCCD_vol_err, "mm^3")
    print(FCCD_vol_err/10**9, "m^3")

    FCCDpctvol_err =FCCD_vol_err/detVol
    print(FCCDpctvol_err*100, " %")



    #calculate FCCD volume - for valentina
    FCCD = 1.9
    #FCCD_vol_1 = np.pi*()*l
    r1, r2 = bottom_rad, bottom_rad-FCCD
    h1, h2, hhat = h, h-FCCD, l
    FCCD_vol_1 = (1/3)*np.pi*(r1**2*h1 - r2**2*h2 - (h1-hhat)*((h1-hhat)*(r1/h1))**2 +(h2-hhat)*((h2-hhat)*(r2/h2))**2)
    FCCD_vol_2 = np.pi*(bottom_rad**2-(bottom_rad-FCCD)**2)*(bottom_height-FCCD) 
    FCCD_vol = FCCD_vol_1 + FCCD_vol_2 + np.pi*((cavity_rad+0.5*FCCD)**2-cavity_rad**2) + np.pi*(cavity_rad+FCCD)**2*(0.5*FCCD)

    FCCDpctvol =FCCD_vol/detVol
    print(FCCDpctvol*100, " %")




    # detector_hits_1 = detector_hits.loc[(detector_hits.z<l+H)&(detector_hits.z>H)]#region 1
    # print("len region 1: ", len(detector_hits_1))
    # detector_hits_FCCD_1 = detector_hits_1.loc[((detector_hits_1.x)**2 + (detector_hits_1.y)**2 < (r**2/h**2)*(detector_hits_1.z-z0)**2)&((detector_hits_1.x)**2 + (detector_hits_1.y)**2 > r_cylinder**2)]
    

    # detector_hits_2 = detector_hits.loc[(detector_hits.z<l2+l+H)&(detector_hits.z>l+H)]#region 2
    # print("len region 2: ", len(detector_hits_2))
    # detector_hits_FCCD_2 = detector_hits_2.loc[((detector_hits_2.x)**2 + (detector_hits_2.y)**2 < B**2)&((detector_hits_2.x)**2 + (detector_hits_2.y)**2 > r_cylinder**2)]

    # detector_hits_3 = detector_hits.loc[(detector_hits.z<l+H+l2+l3)&(detector_hits.z>l+H+l2)]#region 3
    # print("len region 3: ", len(detector_hits_3))
    # detector_hits_FCCD_3 = detector_hits_3.loc[((detector_hits_3.x)**2 + (detector_hits_3.y)**2 < B**2)]

    # detector_hits_FCCD_LIST = [detector_hits_FCCD_1, detector_hits_FCCD_2, detector_hits_FCCD_3]





if __name__ == "__main__":
    main()