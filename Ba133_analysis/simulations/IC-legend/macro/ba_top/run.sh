#!/usr/bin/env bash


cd /lfs/l1/legend/detector_char/enr/hades/simulations/IC-legend/IC160A/Am241/collimated/top/macro

for f in *.mac; do 
	qsub run-g4simple_am_top.qsub  $f;
done 
