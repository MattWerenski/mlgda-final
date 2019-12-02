### This file is meant to serve as an example workflow for analyzing tortuosity of individual vessels.

##  First set your working directory to a location consisting of the .dat files and your directory of R functions

setwd("/Users/alexwork/odrive/Google Drive bruinmail/Research/ucla/savage_lab_vascular_build/tortuosity_calcs_R/code/final/code_and_sample_data")

source("./dat_file_reader.R")
source("./single_vessel_plotter.R")
source("./vessel_poly_fit.R")
source("./vessel_spline_fit.R")
source("./frenet_vectors.R")
source("./tortuosity_metrics.R")
source("./curve_torse_check_plotter.R")
source("./subsample.R")

library("rgl")
library("mgcv")
library("nat")

## Use the dat_file_reader() function to read in a .dat file for analysis
setwd("/Users/alexwork/odrive/Google Drive bruinmail/Research/ucla/savage_lab_vascular_build/extern_materials/audrey_and_kaitlin_2018/angicart_kaitlin/dat files/")

## Vessel ID 50 from next file has a good example of a true torsion spike, and elimination of false torsion spikes.
filename <- "MCAO day 7 microspheres-pecam peri-infarct 1_vd_410x410x1000_561_sm_th_0.37.dat"

## 35 from this one.
filename <- "MCAO day 7 microspheres-pecam NIH 1_vd_410x410x1000_561_sm_th_0.32.dat"

## Use dat_file_reader to convert .dat file to vessels_slice data frame consisting of all vessels with ID listed in order of appearance.
vessels_slice <- dat_file_reader(dat_filename = filename)

## Run the View() function to see the structure of the vessels_slice data.frame.  Note that it contains rescaled coordinates for the vessel backbones (in units of micrometers) and the fourth column is a factor for vessel ID numbers.  Note that these IDs are NOT the same as the vessel nodeids from your other .tsv files.  Integration of the correct vessel ID names happens after performing tortuosity analyses.

View(vessels_slice)

## Quickly use the RGL packages plot3d function to get a coarse view of all vessels in the slice.  The nopen3d() function from the "nat" package allows for panning the 3d plot by pressing down on the trackpad with two fingers and swiping.
nopen3d()
plot3d(vessels_slice[,1:3], col = vessels_slice$ID)

## Alternatively, use the custom made single_vessel_plotter() function for viewing individual vessels.  Index via the which function to capture individual vessels based on their ID number.

# Extracting an individual vessel.  Note that we have also indexed so as to remove the ID column from our single vessel.
vessel <- vessels_slice[which(vessels_slice$ID == 50), 1:3]

nopen3d()
single_vessel_plotter(vessel_coords = as.matrix(vessel), centered = TRUE, frenet = FALSE, new = TRUE)

## Run the new version of vessel_spline_fit using the spline parameter "pspline" for the method that automates the task of identifying the ideal fit.  In short, the spline method itself determines the ideal smoothing parameter via a leave-one-out cross-validation method for a given number of spline knots, spline orders, and penalization degree.  Keep the "m" parameter set at m(4,2).  This sets the exponents of the spline polynomials being used, which in turn influences the number of derivatives to hich continuity is maintained.  The "4" means that the derivatives should be continuous up to and through the third.  The "2" in the "m" parameter sets the penalization degree, which we will also keep fixed at 2.  See Eilers and Marx, "Flexible Smoothing with B Splines and Penaltie".  The method we have coded adds onto this a process to increment the number of knots and test the subsequent fits against one another using an AIC score comparison.  This provides a semi-controlled manner of identifying the ideal fit, but we still have to set a tolerance on the fit itself.  This is done using the "aic_slope" value.  We think that a value of 1 is ideal, but we would like to systematically examine this by varying the aic_slope.  Try integer values of 1, 2, 3, 4, 5.  Note that we can also set the number of knots manually should we want.  This is done by specifying "knots = [integer]".

## Next, when calculating tortuosity metrics, we have a method for controlling for torsion spikes.  This method searches for zeros in the quantity |r'(s)xr''(s)| and filters the spikes when found.  It's a slightly clunky method, but it seems sufficient for now.  The current big problem with this approach is that "zero" on computers are never exactly zero, so we have to set a "filter_tolerance" value that identifies which values are near to zero.  This number may vary with the "subsample_density" value that can also be set.  Defaults for these values are "subsample_density = 2" and "filter_tolerance = 0.01".  For shorter vessels, the subsample_density may need to go higher to have a useful number of points, typically 10 is sufficient.  When varying the "filter_tolerance", try doing so by orders of magnitude.  So other options to try would be 0.1, 0.01, and 0.001.  Feel free to try halfway points as well, such as 0.05, and 0.005.

##The bulk of this next project will be examining the tortuosity measurements under the following conditions.  (i) Experiment with varying "aic_slope", and "filter_tolerance", and as needed "subsample_density", only for the vessels that were previously studied during the summer, so that we can compare the pspline method to the previous cases of identifying smoothing parameters by hand.  (ii) After completing part (i), we should schedule to meet and discuss the results.  But, in the meantime, pick the parameters that seemed to perform best, and experiment with the "scale_up_tortuosity_measures_whole_slice_worksheet.R" file to bulk/batch analyze whole slices for the NIH, IH, and PI regions.  By experiment, I literally mean get creative with how to analyze and compare the measures of tortuosity in the NIH, IH, and PI groups.  Up until now we have mainly been testing our methods to improve them using small samples from these three groups (and there is plenty more testing and improving of our methods to do).  However, I think we are at a good place to start using our methods to look for large scale differences between the groups.  Some suggested plots to consider, (iia) histograms of different tortuosity metrics, (iib) xy-scatter plots of different tortuosity metrics (for example, MT vs. MC) to look for correlation/covariation, (iic) coloring all vessels in a slice by the curvature and torsion metrics.


### Following are some exampes of using the newer functionality.

# Here we try the pspline method of interpolation with m = c(4,2), and no sub-sampling.
vessel_spline_fit(vessel = vessel, number_samples = 10000, spline = "pspline", m = c(4,2), plot = TRUE)

# Note that we have also added versatility for selection of a subsample_density.  For now, we will use subsample_density = 2, which should approximately mean 2 points per micrometer.  
smth_vessel <- vessel_spline_fit(vessel = vessel, number_samples = 100000, spline = "pspline", m = c(4,2), plot = TRUE, subsample_density = 2)

## Now that we have smoothed interpolated vessels, we can add the Frenet vectors to the visualizaitons if we like when using single_vessel_plotter().

nopen3d()
single_vessel_plotter(vessel_coords = smth_vessel[[4]], centered = TRUE, frenet = FALSE, scale = 1.0, col = "purple", frenet_vectors = seq(from = 1, to = nsrow(smth_vessel[[4]]), by = 1))

## Next we'll make graphs of curvature, torsion, and error in curvature and torsion as functions of normalized vessel arc length, using the curve_torse_check_plotter() function.

# To graph curvature and torsion, use plot_type = 1
# Here we experiment using the filter_tolerance.  Note that filtering can be toggled on and off using the command "filter_torsion_spikes = True/False".
curve_torse_check_plotter(vessel_coords = smth_vessel, plot_type = 1, type = "l", filter_torsion_spikes = FALSE, filter_tolerance = 0.01)
curve_torse_check_plotter(vessel_coords = smth_vessel, plot_type = 1, type = "l", filter_torsion_spikes = TRUE, filter_tolerance = 0.01)

# To graph curvature and error in curvature, use plot_type = 2.  Note that these "errors" are not reflective of locations of torsion spikes.
curve_torse_check_plotter(vessel_coords = smth_vessel, plot_type = 2, filter_torsion_spikes = TRUE)

# To graph torsion and error in torsion, use plot_type = 3.  Note that these "errors" are not reflective of locations of torsion spikes.
curve_torse_check_plotter(vessel_coords = smth_vessel, plot_type = 3, filter_torsion_spikes = TRUE)

## To extract various measures of tortuosity, we can use a combination of the functions curvature_torsion_calculator(), distance_metric(), inflection_count_metric, and sum_of_all_angles_metric().  The former of these four methods measures produces multiple metrics all directly related to curvature and torsion.  These are: total and average and maximum curvature and torsion, as well as total and average combined curvature and torsion.  These metrics are discussed in O'Flynn et al., 2007.  The latter three metrics of distance, inflection count, and sum of all angles are discussed in Bullitt et al., 2003.

cur_tor_met <- curvature_torsion_calculator(tangent = smth_vessel[[1]], normal = smth_vessel[[2]], binormal = smth_vessel[[3]], vessel_coords = smth_vessel[[4]], filter_torsion_spikes = TRUE, filter_tolerance = 0.01)

cur_tor_met$curvature
cur_tor_met$torsion
cur_tor_met$TC
cur_tor_met$AC
cur_tor_met$TT
cur_tor_met$AT
cur_tor_met$MC
cur_tor_met$MT
cur_tor_met$TCC
cur_tor_met$ACC

## Distance metric returns both the arclength and the arclength divided by the end-to-end distance
dist_met <- distance_metric(vessel_coords = smth_vessel[[4]])

## Inflection count metric can be toggled to return either the number of inflection points, or the values of the (Delta N)**2 term as a function of vessel arclenth which indicate inflection points when (Delta N)**2 > 2.  Note that this method is likely mis-identifying true inflection points as the above formula is not well defined for 3D curves, but only 2D curves.
inflec_met <- inflection_count_metric(vessel_coords = smth_vessel[[4]], return_count = TRUE, return_coords = FALSE)
inflec_met <- inflection_count_metric(vessel_coords = smth_vessel[[4]], return_count = FALSE, return_coords = TRUE)

## Sum of all angles metric returns the sum of all angles as defined in Bullitt et al., 2003.  Note that this method, adopted by many, is incorrectly calculating spikes in the torsion values if not handled appropriately.

soam <- sum_of_all_angles_metric(smth_vessel[[4]])

## Finally, we can plot some tortuosity metrics on the vessels them selves using color coding.  Here are two examples using curvature and torsion to color the vessel.


# Coloring by curvature.
nopen3d()
single_vessel_plotter(vessel_coords = smth_vessel[[4]], centered = TRUE, frenet = FALSE, col_metric = cur_tor_met$curvature)


# Coloring by torsion
nopen3d()
single_vessel_plotter(vessel_coords = smth_vessel[[4]], centered = TRUE, frenet = TRUE, col_metric = cur_tor_met$torsion, scale = 0.05, frenet_vectors = seq(from = 1050, to = 1150, by = 10))



## Use the above code to experiment with different different fitting parameters for different vessels (and from different slices).  More detailed instruction will come after the weekend.
