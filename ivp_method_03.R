## Separate attempt to use backward Euler with Newton's method for FS equations, then RK4 for position vector since it appeared to work before


library(limSolve)
library(deSolve)
library(rgl)
library(nat)


ivp_method <- function(smth_vessel, arclength = NULL, curvature = NULL, torsion = NULL, tolerance1 = 0.01, mkplot = FALSE, projection = TRUE, seed = 1234, zoom = 1, userMatrix = diag(4), windowRect = c(100, 100, 612, 612), filter_torsion_spikes = FALSE, filter_tolerance = 0.00001, axes = FALSE, aspect = "iso", k = NULL, ...){
  
  
  if(length(smth_vessel) != 4){
    smth_vessel <- frenet_frame_calc(vessel_coords = smth_vessel)
  }
  
  if(is.null(arclength)){
    cur_tor_met_output <- curvature_torsion_calculator(vessel_coords = smth_vessel[[4]], filter_torsion_spikes = filter_torsion_spikes, filter_tolerance = filter_tolerance, k = k)
    arclength <- cur_tor_met_output$arclength
    curvature <- cur_tor_met_output$curvature
    torsion <- cur_tor_met_output$torsion
  }
  
  frntsrt <- function(y, kappa, tau){
    with(as.list(y), {
      dt1 <- kappa*y[4]
      dt2 <- kappa*y[5]
      dt3 <- kappa*y[6]
      
      dn1 <- -kappa*y[1] + tau*y[7]
      dn2 <- -kappa*y[2] + tau*y[8]
      dn3 <- -kappa*y[3] + tau*y[9]
      
      db1 <- -tau*y[4]
      db2 <- -tau*y[5]
      db3 <- -tau*y[6]
      return(c(dt1, dt2, dt3, dn1, dn2, dn3, db1, db2, db3))
    })
  }
  
  jacob <- function(step, kappa, tau){
    bemat <- matrix(data = c(1, 0, 0, -step*kappa, rep(0, 5),
                             0, 1, 0, 0, -step*kappa, rep(0, 4),
                             0, 0, 1, 0, 0, -step*kappa, rep(0, 3),
                             step*kappa, 0, 0, 1, 0, 0, -step*tau, 0, 0,
                             0, step*kappa, 0, 0, 1, 0, 0, -step*tau, 0,
                             0, 0, step*kappa, 0, 0, 1, 0, 0, -step*tau,
                             rep(0, 3), step*tau, 0, 0, 1, 0, 0,
                             rep(0, 4), step*tau, 0, 0, 1, 0,
                             rep(0, 5), step*tau, 0, 0, 1), nrow = 9, ncol = 9, byrow = FALSE)
    return(inv(bemat))
  }
  
  
  ## Remove NaNs from arclength, curvature, and torsion.
  
  nanvec <- which(is.na(curvature))
  arclength <- arclength[-nanvec]
  curvature <- curvature[-nanvec]
  torsion <- torsion[-nanvec]
  smth_vessel$tangent <- smth_vessel$tangent[-nanvec,]
  smth_vessel$normal <- smth_vessel$normal[-nanvec,]
  smth_vessel$binormal <- smth_vessel$binormal[-nanvec,]
  smth_vessel$vessel_coords <- smth_vessel$vessel_coords[-nanvec,]

  ## Initialize output solution for Frenet-Serrat Integration
  phi <- mat.or.vec(nr = length(arclength), nc = 9)
  phi[1,] <- c(smth_vessel$tangent[1,], smth_vessel$normal[1,], smth_vessel$binormal[1,])
  
  for(i in 2:nrow(phi)){
    # print(i)
    step <- arclength[i] - arclength[i-1]
    phi0 <- phi[i-1,]
    jacobinv <- jacob(step = step, kappa = curvature[i], tau = torsion[i])
    phi1 <- phi0 - jacobinv%*%(phi0 - step*frntsrt(y = phi0, kappa = curvature[i], tau = torsion[i])-phi[i-1,])
    print(i)
    while(norm(x = phi1 - phi0, type = "2") > tolerance1){
      norm0 <- norm(x = phi1 - phi0, type = "2")
      phi0 <- phi1
      phi1 <- phi0 - jacobinv%*%(phi0 - step*frntsrt(y = phi0, kappa = curvature[i], tau = torsion[i])-phi[i-1,])
      norm1 <- norm(x = phi1 - phi0, type = "2")
      if(norm1 > norm0){
        print("diverging")
        break
      }
    }
    phi[i,] <- phi1
  }
  
  ## Initialize output solution for position vector.
  pos <- mat.or.vec(nr = length(arclength), nc = 3)
  
  ## Input initial condition for position vector.
  pos[1,] <- smth_vessel$vssl_coords[1,]
  
  ## Execute for loop using Forward Euler
  for(i in 2:nrow(pos)){
    step <- arclength[i] - arclength[i-1]
    pos[i,] <- pos[i-1,] +  step*phi[i-1, 1:3]
  }
  
  out <- list(pos[,1], pos[,2], pos[,3], phi[,1], phi[,2], phi[,3], phi[,4], phi[,5], phi[,6], phi[,7], phi[,8], phi[,9])
  
  # Plot reconstructed vessel (in blue) and original vessel (in black)
  if(mkplot){
    if(projection){
      nopen3d(zoom = zoom, userMatrix = userMatrix, windowRect=windowRect)
      
      plot3d(x = out[[1]], y = out[[2]], z = out[[3]], type = "p", col = "blue", xlab = "", ylab = "", zlab = "", axes = axes, ...)
      plot3d(smth_vessel$vssl_coords, type = "p", col = "black", add = TRUE, ...)
      
      aspect3d(aspect)
      
      box3d()
      
      # Add x grid and projection
      grx <- grid3d('x')
      plot3d(cbind(rgl.attrib(grx[1],'vertices')[1,1], smth_vessel$vssl_coords[,c("y", "z")]),col='gray',add=T, size = 1)
      plot3d(cbind(rgl.attrib(grx[1],'vertices')[1,1], out[[2]], out[[3]]),col='lightskyblue3',add=T, size = 1)
      
      ## Add y grid and projection
      gry <- grid3d('y')
      plot3d(cbind(smth_vessel$vssl_coords[,c("x")],rgl.attrib(gry[1],'vertices')[1,2], smth_vessel$vssl_coords[,c("z")]),col='gray',add=T, size = 1)
      plot3d(cbind(out[[1]], rgl.attrib(gry[1],'vertices')[1,2], out[[3]]),col='lightskyblue2',add=T, size = 1)
      
      ## Add z grid and projection
      grz <- grid3d('z')
      plot3d(cbind(smth_vessel$vssl_coords[,c("x", "y")],rgl.attrib(grz[1],'vertices')[1,3]),col='gray',add=T, size = 1)
      plot3d(cbind(out[[1]], out[[2]],rgl.attrib(grz[1],'vertices')[1,3]),col='lightskyblue3',add=T, size = 1)
      
    }else{
      nopen3d(zoom = zoom, userMatrix = userMatrix, windowRect=windowRect)
      plot3d(x = out[[1]], y = out[[2]], z = out[[3]], type = "p", col = "blue", xlab = "", ylab = "", zlab = "", axes = axes, ...)
      plot3d(smth_vessel$vssl_coords, type = "p", col = "black", add = TRUE, ...)
      box3d()
      grx <- grid3d('x')
      gry <- grid3d('y')
      grz <- grid3d('z')
    }
  }
  
  # Calculate squared error between real and reconstructed vessels.
  magdiff <- mat.or.vec(length(out[[1]]), 1)
  for(i in 1:length(out[[1]])){
    xdiff <- abs(out[[1]][i] - smth_vessel$vssl_coords[i, 1])
    ydiff <- abs(out[[2]][i] - smth_vessel$vssl_coords[i, 2])
    zdiff <- abs(out[[3]][i] - smth_vessel$vssl_coords[i, 3])
    magdiff[i] <- sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
  }
  
  return(list(out, magdiff, sum(magdiff)))
  
}
