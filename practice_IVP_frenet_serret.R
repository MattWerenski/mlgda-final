### Here we experiment with deSolve package for recreating position vector of vessel given initial values and curvature and torsion.

## Mirroring apprach taken in https://stackoverflow.com/questions/42743536/rsolving-ode-with-desolve-package-using-matrices-as-input where we use the ode() function and use matrix inputs to avoid workin with ode.3d(), and the approach taken in https://www.researchgate.net/post/How_can_I_solve_a_system_of_ODEs_with_time_dependent_parameters_in_R_or_in_Python to account for time-varying parameters.

library(deSolve)
library(rgl)
library("nat")

# eqRG = function(tm, state, parms)
# {
#   with(as.list(c(tm, state, parms)),
#        {
#          a1 = parms[["a1"]]
#          a2 = parms[["a2"]]
#          a3 = parms[["a3"]]
#          b1 = parms[["b1"]]
#          b2 = parms[["b2"]]
#          b3 = parms[["b3"]]
#          c1 = parms[["c1"]]
#          c2 = parms[["c2"]]
#          c3 = parms[["c3"]]
#          dy1 = a1(tm) * y2 - b1(tm) * y3 - c1(tm) * y1
#          dy2 = a2(tm) * y1 - b2(tm) * y1 - c2(tm) * y3
#          dy3 = a3(tm) * y2 - b3(tm) * y3 - c3(tm) * y1
#          return(list(c(dy1, dy2, dy3)))
#        })
# }
# tm = seq(0, 10, len = 100)
# state = c(y1 = 1, y2 = 0.5, y3 = 0.02)
# a1 = b1 = c1 = approxfun(tm, - tm/10)
# a2 = b2 = c2 = approxfun(tm, tm * 2)
# a3 = b3 = c3 = approxfun(tm, sin(tm / 20))
# P = list(a1 = a1, a2 = a2, a3 = a3,
#          b1 = a1, b2 = b2, b3 = b3,
#          c1 = c1, c2 = c2, c3 = c3)
# sol = ode(y = state, times = tm, parms = P, func = eqRG)
# plot(sol)


MODEL <- function(s, state, parameters){
  with(as.list(c(state, parameters)), {
    
    kappa <- parameters[["kappa"]]
    tau <- parameters[["tau"]]
    
    r1 <- state["r1"]
    r2 <- state["r2"]
    r3 <- state["r3"]
    t1 <- state["t1"]
    t2 <- state["t2"]
    t3 <- state["t3"]
    n1 <- state["n1"]
    n2 <- state["n2"]
    n3 <- state["n3"]
    b1 <- state["b1"]
    b2 <- state["b2"]
    b3 <- state["b3"]
    
    dr1 <- t1
    dr2 <- t2
    dr3 <- t3
    
    dt1 <- kappa*n1
    dt2 <- kappa*n2
    dt3 <- kappa*n3
    
    dn1 <- -kappa*t1 + tau*b1
    dn2 <- -kappa*t2 + tau*b2
    dn3 <- -kappa*t3 + tau*b3
    
    db1 <- -tau*n1
    db2 <- -tau*n2
    db3 <- -tau*n3
    
    return(list(c(dr1, dr2, dr3, dt1, dt2, dt3, dn1, dn2, dn3, db1, db2, db3)))
  })
}

MODEL(s, state, parameters)

s <- seq(0, 100, length.out = 10000)
# kappa <- c(rep(0.5, 100))
# tau <- rep(1, 100)
kappa <- 10
tau <- 2
parameters <- list(kappa = kappa, tau = tau)
state <- c(r1 = 0, r2 = 0, r3 = 0, t1 = 1, t2 = 0, t3 = 0, n1 = 0, n2 = 1, n3 = 0, b1 = 0, b2 = 0, b3 = 1)

out <- as.data.frame(ode(y = state, times = s, func = MODEL, parms = parameters, "rk4"))
plot3d(x = out$r1, y = out$r2, z = out$r3, type = "l", ylim = c(0, 0.2))


### Here we experiment with writing our own RK4 method so we can incorporate time-dependent parameters.

## myrk4: each element of solnList is a vector
## which contains the grid values
## of one of the dependent variables.
myrk4 = function(init, grid, func) {
  num.var = length(init)
  solnList = list()
  n.grid = length(grid)
  for (k in 1:num.var) {
    solnList[[k]] = vector(length = n.grid)
    solnList[[k]][1] = init[k] }
  h = grid[2] - grid[1] # step size
  yL = vector(length = num.var) # solution at beginning of each step
  for (j in 1:(n.grid-1)) {
    for (k in 1:num.var) yL[k] = solnList[[k]][j]
    k1 = func(grid[j], yL, j) # vector of derivatives
    k2 = func(grid[j] + h/2, yL + h*k1/2, j) # vector of derivatives
    k3 = func(grid[j] + h/2, yL + h*k2/2, j) # vector of derivatives
    k4 = func(grid[j] + h, yL + h*k3, j) # vector of derivatives
    for (k in 1:num.var) solnList[[k]][j+1] = solnList[[k]][j] +
      h*(k1[k] + 2*k2[k] + 2*k3[k] + k4[k])/6 }
  if (num.var==1) solnList[[1]] else solnList
}

frntsrt <- function(t, y, indx){
  with(as.list(y), {
    dr1 <- y[4]
    dr2 <- y[5]
    dr3 <- y[6]
    
    dt1 <- kappa[indx]*y[7]
    dt2 <- kappa[indx]*y[8]
    dt3 <- kappa[indx]*y[9]
    
    dn1 <- -kappa[indx]*y[4] + tau[indx]*y[10]
    dn2 <- -kappa[indx]*y[5] + tau[indx]*y[11]
    dn3 <- -kappa[indx]*y[6] + tau[indx]*y[12]
    
    db1 <- -tau[indx]*y[7]
    db2 <- -tau[indx]*y[8]
    db3 <- -tau[indx]*y[9]
    return(c(dr1, dr2, dr3, dt1, dt2, dt3, dn1, dn2, dn3, db1, db2, db3))
  })
}

tL <- seq(0, 100, length.out = 1000)
kappa <- rep(0.5, 1000)
tau <- rep(0.5, 1000)
tau <- c(seq(0.5, -0.5, length.out = 1000))
tau <- c(seq(0.5, 0.5, length.out = 480), seq(0.5, 10, length.out = 20), seq(10, 0.5, length.out = 20), seq(0.5, 0.5, length.out = 480))
yini <- c(r1 = 0, r2 = 0, r3 = 0, t1 = 1/sqrt(2), t2 = 1/sqrt(2), t3 = 0, n1 = -1/sqrt(2), n2 = 1/sqrt(2), n3 = 0, b1 = 0, b2 = 0, b3 = 1)

out <- myrk4(init = yini, grid = tL, func = frntsrt)
nopen3d()
plot3d(x = out[[1]], y = out[[2]], z = out[[3]], type = "l", col = c(rep("blue", 480), rep("green", 40), rep("red", 500)))

plot(x = 1:1000, y = tau, type = "l")




### Quick practice with real vessel information.


yini <- c(smth_vessel$vssl_coords[5,], smth_vessel$tangent[5,], smth_vessel$normal[5,], smth_vessel$binormal[5,])
tL <- cur_tor_met$arclength[-c(1:4, (nrow(smth_vessel$vssl_coords)-4):nrow(smth_vessel$vssl_coords))]
kappa <- cur_tor_met$curvature[-c(1:4, (nrow(smth_vessel$vssl_coords)-4):nrow(smth_vessel$vssl_coords))]
tau_no_filter <- cur_tor_met$torsion[-c(1:4, (nrow(smth_vessel$vssl_coords)-4):nrow(smth_vessel$vssl_coords))]
tau_with_filter <- cur_tor_met$torsion[-c(1:4, (nrow(smth_vessel$vssl_coords)-4):nrow(smth_vessel$vssl_coords))]
tau <- tau_with_filter

out <- myrk4(init = yini, grid = tL, func = frntsrt)

nopen3d()
plot3d(x = out[[1]], y = out[[2]], z = out[[3]], type = "l", col = "blue")
plot3d(smth_vessel$vssl_coords, type = "l", col = "black", add = TRUE)

nopen3d()
single_vessel_plotter(vessel_coords = as.matrix(data.frame(x = out[[1]], y = out[[2]], z = out[[3]])), centered = FALSE, frenet = FALSE, scale = 1.0, col = "purple", frenet_vectors = seq(from = 1, to = nsrow(smth_vessel[[4]]), by = 1))
single_vessel_plotter(vessel_coords = as.matrix(data.frame(x = out[[1]], y = out[[2]], z = out[[3]])), centered = FALSE, frenet = FALSE, scale = 1.0, col = "purple", frenet_vectors = seq(from = 1, to = nsrow(smth_vessel[[4]]), by = 1))
nopen3d()
single_vessel_plotter(vessel_coords = smth_vessel$vssl_coords, centered = TRUE, frenet = FALSE, scale = 1.0, col = "black", frenet_vectors = seq(from = 1, to = nsrow(smth_vessel[[4]]), by = 1))


magdiff <- mat.or.vec(length(out[[1]]), 1)
for(i in 1:length(out[[1]])){
  xdiff <- abs(out[[1]][i] - smth_vessel$vssl_coords[i+4, 1])
  ydiff <- abs(out[[2]][i] - smth_vessel$vssl_coords[i+4, 2])
  zdiff <- abs(out[[3]][i] - smth_vessel$vssl_coords[i+4, 3])
  magdiff[i] <- sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
}
plot(cur_tor_met$arclength[-c(1:4, (nrow(smth_vessel$vssl_coords)-4):nrow(smth_vessel$vssl_coords))], magdiff, type = "l")
head(magdiff)
max(magdiff)


#### Here we experiment with implementing code to perform error checking within the RK4 IVP solver and fix the torsional values.

helixvssl1D_jit <- as.data.frame(helixvssl1D_jit)

names(helixvssl1D_jit) <- c("x", "y", "z")

smth_vessel <- vessel_spline_fit(vessel = helixvssl1D_jit, number_samples = 100000, spline = "pspline", aic_slope = 5, m = c(4,2), plot = TRUE, subsample_density = 100)

curve_torse_check_plotter(vessel_coords = smth_vessel, plot_type = 1, type = "l", filter_torsion_spikes = FALSE, filter_tolerance = 0.01, save_plot = FALSE)
cur_tor_met <- curvature_torsion_calculator(tangent = smth_vessel[[1]], normal = smth_vessel[[2]], binormal = smth_vessel[[3]], vessel_coords = smth_vessel[[4]], filter_torsion_spikes = FALSE, filter_tolerance = 0.05)


yini <- c(smth_vessel$vssl_coords[5,], smth_vessel$tangent[5,], smth_vessel$normal[5,], smth_vessel$binormal[5,])
tL <- cur_tor_met$arclength[-c(1:4, (nrow(smth_vessel$vssl_coords)-4):nrow(smth_vessel$vssl_coords))]
kappa <- cur_tor_met$curvature[-c(1:4, (nrow(smth_vessel$vssl_coords)-4):nrow(smth_vessel$vssl_coords))]
tau <- cur_tor_met$torsion[-c(1:4, (nrow(smth_vessel$vssl_coords)-4):nrow(smth_vessel$vssl_coords))]

tau[4501:4550] <- seq(from = tau[4501], to = -10, length.out = 50)
tau[4551:4600] <- seq(from = -10, to = tau[4600], length.out = 50)

out <- myrk4(init = yini, grid = tL, func = frntsrt)


nopen3d()
single_vessel_plotter(vessel_coords = as.matrix(data.frame(x = out[[1]], y = out[[2]], z = out[[3]])), centered = TRUE, frenet = FALSE, scale = 1.0, col = "purple", frenet_vectors = seq(from = 1, to = nsrow(smth_vessel[[4]]), by = 1))


magdiff <- mat.or.vec(length(out[[1]]), 1)
for(i in 1:length(out[[1]])){
  xdiff <- abs(out[[1]][i] - smth_vessel$vssl_coords[i+4, 1])
  ydiff <- abs(out[[2]][i] - smth_vessel$vssl_coords[i+4, 2])
  zdiff <- abs(out[[3]][i] - smth_vessel$vssl_coords[i+4, 3])
  magdiff[i] <- sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff)
}

plot(cur_tor_met$arclength[-c(1:4, (nrow(smth_vessel$vssl_coords)-4):nrow(smth_vessel$vssl_coords))], tau)
plot(cur_tor_met$arclength[-c(1:4, (nrow(smth_vessel$vssl_coords)-4):nrow(smth_vessel$vssl_coords))], magdiff, type = "l")

thresh <- 0.05

##  Major deviations between curves happen AFTER spikes in torsion occur.  Will need to back check and correct, then recalculate.  Once deviation is detected during main routine, start subroutine and use r'Xr'' to search for most recent minimum in curvature/spike in torsion, then correct value of torsion to zero, and restart solver process, returning to main routine.



