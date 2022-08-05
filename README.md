# Failure_surface_estimator
## Method to estimate the initial landslide failure surface and volumes using grid points and spline curves in MATLAB
Failure_surface_estimator is a new method to estimate the initial landslide failure surface and volumes using grid points and spline curves in MATLAB. The data inputs and the MATLAB codes and functions are discribed in the main code as instructions while running, to complete all the steps.profles. The model will give the depth of the probable failure surface plotted using the 2D grid function.
We can easily visualize the results by running the codes; the display will show the step-by-step processes involved.The user needs to have the MATLAB Mapping Toolbox (MATLAB 2021) installed as an extension to run the code.This is compatible with organizing geographic data and allows for the interpolation, trimming, resampling, and transformation of coordinates.

There are many functions required to compile the codes together to get the results. It  includes -deg2utm, kml2struct, gridft, poly3n, and splin_param_new. For smooth running, better to put all the functions and data at same location in a single file.  

## Data required
The following data are required:

1. KML fle of the contour lines tracing the boundary of the failure surface

2. DEM of the area in TIFF format covering the contour limits

### QUICK START:

Open Matlab, select the Main_code.m script, run the script and visualise the results as matlab figures and read the valuse from  commond window. For visualising the results step by step, select the main_code.mlx script and run it. You will see the outputs on the right side. 

### COMPATIBILITY:

Failure_surface_estimetor is designed under a MATLAB architecture. However, it can also be run with OCTAVE. We did not check the compatibility of the estimetor with every version of MATLAB nor OCTAVE, but we provide a non-exhaustive list of compatibility.

MATLAB version:R2020b,R2020a, R2018b, R2018a, R2017b, R2016a, R2013b

OCTAVE version: 5.1.0.0

## Funding
Open access funding provided by University of Lausanne - Switzerland.

### Contact:gautam@es.iitr.ac.in
