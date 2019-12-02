vessel_poly_fit <- function(vessel, number_samples = 10000, poly_degree = 5, plot = TRUE, subsample_density, ...){
  ## Generic comments about this function...
  
  library("pracma")
  
  vessel_steps <- nrow(vessel)
  
  if(plot){
    # source("./single_vessel_plotter.R")
    nopen3d()
    single_vessel_plotter(vessel_coords = vessel, col = "black", ...)
  }
  
  ## Here begins the polynomial fitting
  ## We need to add a column indexing the voxels.  This appears to be needed to the way that the predict() function identifies and tracks data.  
  fit_length <- number_samples
  
  vessel <- data.frame(vessel, "voxel" = 1:vessel_steps)
  
  fitx <- lm(formula = x ~ poly(voxel, degree = poly_degree, raw = TRUE), data = vessel)
  fity <- lm(formula = y ~ poly(voxel, degree = poly_degree, raw = TRUE), data = vessel)
  fitz <- lm(formula = z ~ poly(voxel, degree = poly_degree, raw = TRUE), data = vessel)
  
  smth_vessel <- data.frame(voxel = seq(from = 1, to = vessel_steps, length.out = fit_length))

  smth_vessel$predx <- predict(fitx, smth_vessel)
  smth_vessel$predy <- predict(fity, smth_vessel)
  smth_vessel$predz <- predict(fitz, smth_vessel)
  
  smth_vessel <- smth_vessel[,-1]
  smth_vessel <- as.matrix(smth_vessel)
  colnames(smth_vessel) <- c("x", "y", "z")
  
  if(!is.null(subsample_density)){
    tic("subsampling")
    subsample_indeces <- subsample(vessel_coords = smth_vessel, density = subsample_density)
    smth_vessel <- smth_vessel[subsample_indeces,]
    toc()
  }

  
  if(plot){
    single_vessel_plotter(vessel_coords = smth_vessel, new = TRUE, col = "blue", ...)
  } 
  
  # Here we initialize the Frenet-Serrat frame vectors.
  tangent_array <- mat.or.vec(nrow(smth_vessel), 3)
  normal_array <- mat.or.vec(nrow(smth_vessel), 3)
  binormal_array <- mat.or.vec(nrow(smth_vessel), 3)
  
  tangent_array[,] <- NaN
  normal_array[,] <- NaN
  binormal_array[,] <- NaN
  print(nrow(smth_vessel))
  for(i in 3:(length(smth_vessel[,1]) - 2)){
    print(i)
    tangent_array[i,] <- tangent_vector(smth_vessel[i-2,], smth_vessel[i-1,], smth_vessel[i,], smth_vessel[i+1,], smth_vessel[i+2,])
    normal_array[i,] <- normal_vector(smth_vessel[i-2,], smth_vessel[i-1,], smth_vessel[i,], smth_vessel[i+1,], smth_vessel[i+2,])
    binormal_array[i,] <- binormal_vector(smth_vessel[i-2,], smth_vessel[i-1,], smth_vessel[i,], smth_vessel[i+1,], smth_vessel[i+2,])
  }
  
  return(list(tangent = tangent_array, normal = normal_array, binormal = binormal_array, vssl_coords = smth_vessel))
}
