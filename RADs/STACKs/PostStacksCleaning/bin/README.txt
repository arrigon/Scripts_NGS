Install:
This pipeline depends on the R packages:
doMC
foreach
multicore

Running: 
The scripts have to process large STACKs outputs. This processing is typically long and lots of time can be gained by re-using intermediate working files. 
These guys are saved in savedanalyses/. Thus leave these files in there if planning to iterate the Postcleaning script several times, with distinct parameters.

