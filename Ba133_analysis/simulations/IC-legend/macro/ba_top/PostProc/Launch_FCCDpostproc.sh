#!/usr/bin/env bash
export PATH=~/miniconda3/bin:$PATH
#MC_file_ID = 'IC160A_ba_top_coll_01'
#python /lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/simulations/IC-legend/macro/ba_top/PostProc/postprochdf5.py detector_IC160A_ba_top_81mmNEW4_01 #IC160A_ba_top_gammas_81mmNEW4_01 #IC160A_ba_top_gammas_81mmNEW3_01 #IC160A_ba_top_gammas_01 #$MC_file_id
python /lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/simulations/IC-legend/macro/ba_top/PostProc/postproc_FCCD.py IC160A_ba_top_81mmNEW4_01 #IC160A_ba_top_gammas_01 #$MC_file_id
#python /lfs/l1/legend/users/aalexander/HADES_detchar/Ba133_analysis/simulations/IC-legend/macro/ba_top/PostProc/drawPostProcessed_FCCD.py IC160A_ba_top_gammas_01 #$MC_file_id


