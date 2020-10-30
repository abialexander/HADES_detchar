#!/bin/bash

PATH_SIMULATION=/lfs/l1/legend/detector_char/enr/hades/simulations/legend-g4simple-simulation/IC-legend/general_geometry
DETECTOR=IC160A
SOURCE=th
POSITION=_top
TABLE=1

ln -sf $PATH_SIMULATION/crystal_$DETECTOR.gdml
#ln -sf $PATH_SIMULATION/teflon_wrap.gdml #old, no longer exists, wrap is not actually made from teflon
ln -sf $PATH_SIMULATION/wrap.gdml
ln -sf $PATH_SIMULATION/holder.gdml
#ln -sf $PATH_SIMULATION/Alcap.gdml
ln -sf $PATH_SIMULATION/cryostat.gdml
ln -sf $PATH_SIMULATION/bottom_plate.gdml
ln -sf $PATH_SIMULATION/lead_castle_table1.gdml
#ln -sf $PATH_SIMULATION/source_encapsulated_$SOURCE.gdml
#ln -sf $PATH_SIMULATION/source_encapsulated_th.gdml #new 12/10/20
ln -sf $PATH_SIMULATION/collimator_Cu$POSITION.gdml

#ln -sf $PATH_SIMULATION/define_source_$SOURCE.xml define_source_th.xml #dont use

ln -sf $PATH_SIMULATION/define_lead_castle.xml

ln -sf $PATH_SIMULATION/define_$DETECTOR.xml define_detector.xml



