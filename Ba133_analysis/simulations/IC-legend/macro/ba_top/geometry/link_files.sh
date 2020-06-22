#!/bin/bash

PATH_SIMULATION=/lfs/l1/legend/detector_char/enr/hades/simulations/legend-g4simple-simulation/IC-legend/general_geometry
DETECTOR=IC160A
SOURCE=am
POSITION=_top
TABLE=1

ln -sf $PATH_SIMULATION/crystal_$DETECTOR.gdml
ln -sf $PATH_SIMULATION/teflon_wrap.gdml
ln -sf $PATH_SIMULATION/holder.gdml
ln -sf $PATH_SIMULATION/Alcap.gdml
ln -sf $PATH_SIMULATION/bottom_plate.gdml
ln -sf $PATH_SIMULATION/lead_castle_table1.gdml
#ln -sf $PATH_SIMULATION/source_encapsulated_$SOURCE.gdml
ln -sf $PATH_SIMULATION/collimator_Cu$POSITION.gdml

ln -sf $PATH_SIMULATION/define_source_$SOURCE$POSITION.xml define_source.xml

ln -sf $PATH_SIMULATION/define_lead_castle.xml

ln -sf $PATH_SIMULATION/define_$DETECTOR.xml define_detector.xml



