 IC_geometry
=============




IC-gerda folder contains the file concerning the geometry and MC simulation of the IC detectors used in Gerda experiment - IC50A, IC50B, IC48A, IC48B, IC74A.

IC-legend folder contains the file concerning the geometry and MC simulation of the other IC detectors  - IC60A, IC60B, IC62B, IC66B.

The dimensions of each detector (Crystal, Holder, Al cap, Positions) are taken from the Mirion documents and are defined in <code>define_IC\*.xml</code>.   
The dimensions of the two tables (Table1, Table2) are defined in <code>define_lead_castle.xml</code>.  
The dimensions of the source are defined in <code>define_source.xml</code>.

The geometry of each component of the detector (crystal, teflon wrap,  holder, Alcap) and of the table (bottom_plate, lead_castle_table1, lead_castle_table2, copper_plate_table2) is defined in each respective GDML module.

<code>detector.gdml</code> is the mother GDML file which includes all the GDML modules.


<code>run_detector_vis.mac</code> is the macro file used for the visualization of the geometry.  
<code>run_detector.mac</code> is the macro file used for statistics.





### GDML files


GDML modules are GDML files used in the definition of geometries. They allow an easier to understand view of big geometries as they split it into smaller (and therefore more readable) pieces.
A GDML module is a normal GDML file and it is completely independent from all the other modules. It contains all the information (definitions, materials, solids and volumes) it needs, to be fully defined.
To include a GDML module in a mother module (<code>detector.gdml</code>), in the structure section of the mother module you should add a physical volume which points to the child module.

The enriched germanium material is manually defined in <code>crystal.gdml</code>; the other materials are taken from Geant4's NIST Material Database (the GDML parser will complain that the materials have not been defined, but Geant4 will still run without error).

The source and the components of the detectors are build using polycones.  
The center of coordinates (0,0,0) of the volume is set at the top of each polycone, while for all other solids (e.g. box,tube,cone) it is set at the center.  
It is important to take in mind that  when positioning a volume associated to a Boolean solid, the relative center of coordinates considered for the positioning is the one related to the first of the two constituent solids.


In order to choose one detector you can use:
<code>ln -sf define_namedetector.xml define_detector.xml</code>  
For IC-legend it is not enough. You have to comment the unused detectors in World structure of <code>detector.gdml</code>.

In IC-legend the distance of the source from the top of Al cap is equal to 4.2 cm for IC160A and 4.5 cm for the other detectors. Then, it is required to set the right value in 
<code>position_source_fromAlcap</code>) in <code>define_source.xml</code> file and in<code>/gps/pos/centre 0 0 -4.2</code> in  <code>run_detector.mac</code> file.  


In World structure of <code>detector.gdml</code> (IC-gerda and IC-legend) you can also choose table1 or table2. 


To draw the GDML file in ROOT:  
<code>TGeoManager::Import("file.gdml");  
gGeoManager->GetTopVolume()->Draw("ogl");</code>



### MACRO files

In MAC file there are all the settings for the simulation (e.g. randomseed set, physics list set, GDML detector set, output set...).  
You can also set the type of the source.

In <code>run_detector_vis.mac</code> the settings for the visualization are included.  
The gamma filter is on.  
Different energy tracks are related to different colors.

The output is set as root file (and heprep file for the visualization).






