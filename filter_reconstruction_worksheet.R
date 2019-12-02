setwd("/Users/alexwork/odrive/Google Drive bruinmail/Research/ucla/savage_lab_vascular_build/tortuosity_calcs_R/code/final/code_and_sample_data")

source("./dat_file_reader.R")
source("./single_vessel_plotter.R")
source("./vessel_poly_fit.R")
source("./vessel_spline_fit.R")
source("./frenet_vectors.R")
source("./tortuosity_metrics.R")
source("./curve_torse_check_plotter.R")
source("./subsample.R")
source("/Users/alexwork/odrive/Google Drive bruinmail/Research/ucla/savage_lab_vascular_build/tortuosity_calcs_R/code/drafts/ivp_method_02_original.R")

library("rgl")
library("mgcv")
library("nat")

setwd("/Users/alexwork/odrive/Google Drive bruinmail/Research/ucla/savage_lab_vascular_build/extern_materials/audrey_and_kaitlin_2018/angicart_kaitlin/dat files/")

## Vessel ID 50 from next file has a good example of a true torsion spike, and elimination of false torsion spikes.  Downside is correct epsilon ball around zero is 10 times bigger than the other examples from different image set.  Still worth trying this approach, coding up somehow for Kaitlin to examine.
filename <- "MCAO day 7 microspheres-pecam peri-infarct 1_vd_410x410x1000_561_sm_th_0.37.dat"

vessels_slice <- dat_file_reader(dat_filename = filename)

## Quickly use the RGL packages plot3d function to get a coarse view of all vessels in the slice.
# nopen3d()
# plot3d(vessels_slice[,1:3], col = vessels_slice$ID)

## Use vessels 12 and 50, 55 is interesting as it is quite straight.
vessel <- vessels_slice[which(vessels_slice$ID == 55), 1:3]

smth_vessel <- vessel_spline_fit(vessel = vessel, number_samples = 20000, spline = "pspline", aic_slope = 5, m = c(4,2), plot = TRUE, subsample_density = 10)

## Plotter for single vessel here.
nopen3d()
single_vessel_plotter(vessel_coords = smth_vessel[[4]], centered = TRUE, frenet = TRUE, scale = 0.75, col = "purple", frenet_vectors = seq(from = 25, to = 35, by = 1))


curve_torse_check_plotter(vessel_coords = smth_vessel[[4]], plot_type = 1, type = "p", filter_torsion_spikes = FALSE, filter_tolerance = 0.05, save_plot = FALSE, main = "Vessel 55, 4.75X Sampling", filename = "vessel_55_4.75X_sampling_cur_tor")

# Use below block of commented code to trim noisy ends of vessels for short-term fix.  Long term fix means add parameter to vessel_spline_fit() that allows for toggling how much of vessels to trim.

start <- 5

smth_vessel <- list(tangent = smth_vessel[[1]][start:nrow(smth_vessel[[1]]),], normal = smth_vessel[[2]][start:nrow(smth_vessel[[1]]),], binormal = smth_vessel[[3]][start:nrow(smth_vessel[[1]]),], vssl_coords = smth_vessel[[4]][start:nrow(smth_vessel[[1]]),])

cur_tor_met <- curvature_torsion_calculator(vessel_coords = smth_vessel[[4]], filter_torsion_spikes = FALSE, filter_tolerance = 0.005)

ivp_output <- ivp_method(smth_vessel = smth_vessel, arclength = cur_tor_met$arclength, curvature = cur_tor_met$curvature, torsion = cur_tor_met$torsion, tolerance1 = 0.1, tolerance2 = 0.1, mkplot = TRUE, projection =TRUE)

rgl.postscript(filename = "/Users/alexwork/Desktop/vessel_55_4.75X_sampling_backeul.pdf", fmt = "pdf")


epar3d(pp)

plot3d()
?par3d()
pp <- par3d(no.readonly = TRUE)
par3d(pp)


#### Sample code for fixing the perspective of plot3d().
# ## In an inital session:
# 
# library(rgl)
# plot3d(iris) 
# 
# ## Now move the image around to an orientation you like
# 
# ## Save RGL parameters to a list object
# pp <- par3d(no.readonly=TRUE)
# 
# ## Save the list to a text file
# dput(pp, file="irisView.R", control = "all")
# 
# .......
# 
# ## Then, in a later session, to recreate the plot just as you had it:
# 
# library(rgl)
# pp <- dget("irisView.R")
# plot3d(iris)
# par3d(pp)
