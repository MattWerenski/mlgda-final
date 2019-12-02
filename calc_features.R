
setwd("D:\\work\\Classes\\Tufts\\MLonGraphs\\Project\\source")
source("tortuosity_metrics.R")
source("vessel_spline_fit.R")



source("dat_file_reader.R") # reads the .dat files output from angicart
source("single_vessel_plotter.R") # specialized plotting function for visualizing individual vessels
source("vessel_poly_fit.R") # function for performing polynomial fits of vessels
source("vessel_spline_fit.R") # function for performing splining fits of vessels
source("frenet_vectors.R") # function for calculating the Frenet-Serret frame vectors
source("tortuosity_metrics.R") # function for calculating many tortuosity metrics
source("curve_torse_check_plotter.R") # function for plotting curvature and torsion versus normalized arclength
source("subsample.R") # function needed for subsampling vessels at the end of the splining procedure.
source("ivp_method_03.R") # function for reconstructing vessels based on measured values of curvature and torsion

library("rgl")
library("mgcv")
#library("nat")

## Use the dat_file_reader() function to read in a .dat file for analysis
# First set the working directory...for example
setwd("D:\\work\\Classes\\Tufts\\MLonGraphs\\Project\\data\\lung01_024")


filename <- "lung01_024_coords_scaled.csv"
## 35 from this one.

##### Warning ##### 
##### The remainder of this tutorial assumes you are working with vessel number 50 from "MCAO day 7 microspheres-pecam peri-infarct 1_vd_410x410x1000_561_sm_th_0.37.dat". ######

# Run the next line of code to use the dat_file_reader() function and define its output as vessels_slice
vessels_slice <- read.csv(filename)

## Run the View() function to see the structure of the vessels_slice data.frame.  Note that it contains rescaled coordinates for the vessel backbones (in units of micrometers) and the fourth column is a factor for vessel ID numbers.  Note that these IDs are NOT the same as the vessel nodeids from your other .tsv files.  Integration of the correct vessel ID names happens after performing tortuosity analyses.

## Alternatively, use the custom made single_vessel_plotter() function for viewing individual vessels.  Index via the which function to capture individual vessels based on their ID number.

# Extracting an individual vessel.  Note that we have also indexed so as to remove the ID column from our single vessel.
# Try  3, 12, 25, 32
#vessel <- vessels_slice[which(vessels_slice$ID == 1), 1:3]
#smth_vessel <- vessel_spline_fit(vessel = vessel, number_samples = 20000, spline = "pspline", m = c(4,2), plot = FALSE, aic = 5, subsample_density = 10)

vessels <- read.csv("D:\\work\\Classes\\Tufts\\MLonGraphs\\Project\\data\\test_lung_withroots\\lung_coords_scaled.csv")
vessels$nodeid<-as.factor(vessels$nodeid)
#levels(vessels)
ids <- unique(vessels$nodeid)


  
feats = c('TC','AC','TT','AT','MC','MT','TCC','ACC', 'TCsq','TTsq','local_minima','comb_curv_int')
feats_name = c('TC','AC','TT','AT','MC','MT','TCC','ACC', 'TCsq','TTsq','local_minima','comb_curv_int','name')
name_feats = c('name','TC','AC','TT','AT','MC','MT','TCC','ACC', 'TCsq','TTsq','local_minima','comb_curv_int')

id = 2



v<-data.matrix(data.frame(vessels$x[vessels$nodeid == ids[id]],vessels$y[vessels$nodeid == ids[id]],vessels$z[vessels$nodeid == ids[id]]))
#v<-vessel_spline_fit(vessel = v, subsample_density = 100, plot = FALSE, m = c(2,2))[[4]]
smth_vessel <- vessel_spline_fit(vessel = v, number_samples = 20000, spline = "pspline", m = c(4,2), plot = FALSE, aic = 5, subsample_density = 10)


metrics <- curvature_torsion_calculator(v)

res<-metrics[feats]
#res<-append(res,as.character(ids[1]) )

#res<-c(ids[1],res)
out<-data.frame(res)
#colnames(out)<-feats_name

setwd("D:\\work\\Classes\\Tufts\\MLonGraphs\\Project\\source")



for (i in 2:length(ids))
{
  
  v<-data.matrix(data.frame(vessels$x[vessels$nodeid == ids[i]],vessels$y[vessels$nodeid == ids[i]],vessels$z[vessels$nodeid == ids[i]]))
  
 
  if(nrow(v) <= 4 )
  {
    #print(ids[i])
    res <-  rep(NaN,length(feats))
  }
  else{
    #print(nrow(v))
    v <- vessel_spline_fit(vessel = v, number_samples = 20000, spline = "smoothing", m = c(4,2), plot = FALSE, aic = 5, subsample_density = 10)
    
    metrics <- curvature_torsion_calculator(v)
    
    res<-metrics[feats]
    print(res)
    #res<-append(res,as.character(ids[i]) )
    #colnames(res)<-feats_name
    out<-rbind(out,res)
  }
  
  
  
}
#out<-out[name_feats]