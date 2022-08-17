~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Adapted Nature-based solutions Rational Method (ANaRM) 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

~~~~~~ Introduction ~~~~~~
R-code plus example data to implement the "ANaRM" method described in Miller et al. 2022. Essentially, the code implements the following equation for peak flow:

   Qp =  Cr * Cv * i * Area * FARL^theta
 
Where "Qp" is the peak flow, "Cr" is a dimensionless routing coefficient set as 1.3 by default, "Cv" is the mean runoff ratio over the upstream "Area", "i" is the mean hourly rainfall over the design storm event duration, and "FARL" (Flood Attenuation by Reservoirs and Lakes) with a default power of theta = 3.445 accounts for the reduction in peak flow due to water bodies. This equation is implemented on a grid to produce a peak flow raster and other useful rasters. It is best used for small urban catchments and especially can be used to compare "GREEN" scenarios (with added nature-based solutions) to a "BASE" scenario (Miller 2022).


~~~~~~ Files ~~~~~~
* ANaRM_main.R : Main R-code.
* ANaRM_functions.R : Functions called by ANaRM_main.R
* options.inp : An example of a basic options file.
* runoff.csv : A table of Cv values for each land class.
* OUTLETS : An example of how locations of interest (flow checkpoints inside the catchment) can be specified.
* LC_28039.tif : A modified version of the 10m resolution open access (CC BY 4.0, https://creativecommons.org/licenses/by/4.0/) ESA Worldcover (Zanaga et al., 2020, https://esa-worldcover.org/en/data-access) for the Rea catchment, Birmingham. See Miller et al. 2022 for details of modifications.
* LIDAR_outf1.tif : Outflow directions, derived from LIDAR Composite DTM 2019 (Open Government Licence v3.0, https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/) using the "Fill" and "Flow Direction" functions in ArcMap. 
!!BE AWARE!! - part of the catchment was not captured in the LIDAR elevation data. These flow directions are provided as an example only. They were not used in Miller et al. 2022, which instead used NEXTMap data (Intermap Technologies, 2007) for the entire study area. 
* LC_GREEN_28039.tif : An example of a "GREEN" scenario based on LC_28039.tif.
* PONDS.tif : An example of how additional ponds may be specified (for a "GREEN" scenario) if useLcmScnForFarlScn = F

~~~~~~ Example usage ~~~~~~

Basic example:
This is a basic example that uses the "BASE" landcover map and the flow directions to create the following:
- flow_base_scenario.tif (peak flows for BASE scenario) in m^3/s.
- rea_checkpoints_example_OUTPUT.csv (flow data at locations specified in the file pointed to by outletsInFN).
- farl_base.tif (gridded FARL values - using the slightly modified definition specified in Miller et al. 2022).
- CCAR.tif (Cumulative Catchment Area) in units of the grid cell resoution (10x10m).

1. Place ANaRM_main.R, ANaRM_functions.R, options.inp and runoff.csv in a base directory. Create subdirectories "input", "output1", and "output2".
2. Place LC_28039.tif, LIDAR_outf1.tif and **OUTLETS** in the "input" directory (place LC_GREEN_28039.tif & PONDS.tif here also, if "GREEN" scenarios are to be tested).
3. Check that the base directory created in step 1 is the working directory (getwd() in R) that R is using, or else edit baseDir in options.inp and optionsFN in ANaRM_main.R. 
4. Run ANaRM_main.R - outputs will be created in the output1 directory.

Advanced usage:
ANaRM is designed to be run in two stages to test multiple "GREEN" scenarios quickly. See Miller et al. 2022 and the comments in options.inp and ANaRM_main.R. LC_GREEN_28039.tif and/or PONDS.tif can be used to test this.   
TO DO: add example two stage options files.

~~~~~~ Questions/comments ~~~~~~
johwal@ceh.ac.uk or wallbank.j.r@gmail.com


~~~~~~References ~~~~~~
Miller, J.D., Vesuviano, G., Wallbank, J.R., Fletcher, D.H. & Jones, L. (2022). Hydrological assessment of urban NBS in Ecosystem Service toolkit applications. Submitted 2022.

Zanaga, D., Van De Kerchove, R., De Keersmaecker, W., et al. (2021). ESA WorldCover 10 m 2020 v100. https://doi.org/10.5281/zenodo.5571936 

LIDAR Composite DTM 2019 10m https://www.data.gov.uk/dataset/8311f42d-bddd-4cd4-98a3-e543de5be4cb/lidar-composite-dtm-2019-10m 

Intermap Technologies (2007): NEXTMap British Digital Terrain Model Dataset Produced by Intermap. NERC Earth Observation Data Centre (Accessed 20/06/22) https://catalogue.ceda.ac.uk/uuid/8f6e1598372c058f07b0aeac2442366d
