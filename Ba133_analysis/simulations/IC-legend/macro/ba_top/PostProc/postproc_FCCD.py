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


    # Modify this value for different energy resolution
    pctResAt1MeV = 0.15 #0.15 #5#;

    hdf5_path = "/lfs/l1/legend/users/aalexander/hdf5_output/"

    if(len(sys.argv) != 2):
        print('Usage: postprochdf5.py [IC160A_ba_top_coll_01]')
        sys.exit()

    MC_file_id = sys.argv[1]
    #MC_file = 'IC160A_ba_top_coll_01'

    # have to open the input file with h5py (g4 doesn't write pandas-ready hdf5)
    g4sfile = h5py.File(hdf5_path+'detector_'+MC_file_id+'.hdf5', 'r')
    g4sntuple = g4sfile['default_ntuples']['g4sntuple']

    g4sdf = pd.DataFrame(np.array(g4sntuple), columns=['event'])

    # # build a pandas DataFrame from the hdf5 datasets we will use
    g4sdf = pd.DataFrame(np.array(g4sntuple['event']['pages']), columns=['event'])
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['step']['pages']), columns=['step']),lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['Edep']['pages']), columns=['Edep']),lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['volID']['pages']),columns=['volID']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['iRep']['pages']),columns=['iRep']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['x']['pages']),columns=['x']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['y']['pages']),columns=['y']), lsuffix = '_caller', rsuffix = '_other')
    g4sdf = g4sdf.join(pd.DataFrame(np.array(g4sntuple['z']['pages']),columns=['z']), lsuffix = '_caller', rsuffix = '_other')

    # apply E cut / detID cut and sum Edeps for each event using loc, groupby, and sum
    # write directly into output dataframe
    detector_hits = g4sdf.loc[(g4sdf.Edep>0)&(g4sdf.volID==1)]
    #print(detector_hits)
    print(detector_hits.size)


    #first save original file, i.e. no FCCD
    # procdf = pd.DataFrame(detector_hits.groupby(['event','volID','iRep'], as_index=False)['Edep'].sum())
    # procdf = procdf.rename(columns={'iRep':'detID', 'Edep':'energy'})
    # # apply energy resolution function
    # procdf['energy'] = procdf['energy'] + np.sqrt(procdf['energy'])*pctResAt1MeV/100.*np.random.randn(len(procdf['energy']))
    # print(procdf.size)
    # procdf.to_hdf(hdf5_path+'processed/processed_detector_'+MC_file_id+'.hdf5', key='procdf', mode='w')

    #apply FCCD cuts
    top_width = 75.5 #mm - these come from detector.xml file
    bottom_width = 79.8 #mm
    #FCCD_list = [0.25]
    FCCD_list = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5] #mm

    for FCCD in FCCD_list:

        AV_radius = top_width/2 - FCCD
        detector_hits_FCCD = detector_hits.loc[((g4sdf.x)**2 + (g4sdf.y)**2 < AV_radius**2)]
        #print(detector_hits_FCCD)
        print(str(FCCD))
        print(detector_hits_FCCD.size)

        procdf = pd.DataFrame(detector_hits_FCCD.groupby(['event','volID','iRep'], as_index=False)['Edep'].sum())
        procdf = procdf.rename(columns={'iRep':'detID', 'Edep':'energy'})
        # apply energy resolution function
        procdf['energy'] = procdf['energy'] + np.sqrt(procdf['energy'])*pctResAt1MeV/100.*np.random.randn(len(procdf['energy']))
        print(procdf.size)
        procdf.to_hdf(hdf5_path+'processed/processed_FCCD'+str(FCCD)+'mm_detector_'+MC_file_id+'.hdf5', key='procdf', mode='w')

if __name__ == "__main__":
    main()