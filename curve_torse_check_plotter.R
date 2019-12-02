curve_torse_check_plotter <- function(vessel_coords, plot_type = 1, type = "p", filter_torsion_spikes = FALSE, filter_tolerance = 0.00001, save_plot = FALSE, format = "png", filename = "curve_torse_vs_arclength", main = NULL, ylimcurv = NULL, ylimtors = NULL, ylabcurv = "Curvature", ylabtors = "Torsion", pointsize = 10, width = 8, height = 6, las = 0, tortuosity_metrics = NULL, k = NULL, ...){
  ## This function serves as an example for plotting curvature/torsion versus normalized arclenth, curvature/curvature-error versus normalized arclength, and torsion/torsion-error versus normalized arclength.
  
  # source("./tortuosity_metrics.R")
  # source("./frenet_vectors.R")
  # 
  if(is.null(tortuosity_metrics)){
    tortuosity_metrics <- curvature_torsion_calculator(vessel_coords = vessel_coords, filter_torsion_spikes = filter_torsion_spikes, filter_tolerance = filter_tolerance, k = k)
  }else{
    tortuosity_metrics <- tortuosity_metrics
    }
  
  

  ## Plotting curvature and torsion in one graph
  if(plot_type == 1){
    if(is.null(main)){
      main = "Curvature and Torsion vs Normalized Arclength"
    }
    if(save_plot){
      if(format == "png"){
        png(filename = paste(filename, ".png", sep = ""), width = width, height = height, res = 300, units = "in")
        par(mfrow = c(2, 1),     # 2x2 layout
            oma = c(5, 1, 1, 1), # 4 rows of text at the outer left and bottom margin
            mar = c(1, 4.1, 2, 0), # space for one row of text at ticks and to separate plots
            mgp = c(2.75, 1, 0), 
            las = las, 
            cex = 1.2) 
        plot(x = tortuosity_metrics$arclength/max(tortuosity_metrics$arclength, na.rm = T), tortuosity_metrics$curvature, xlab = NULL, ylab = ylabcurv, axes = FALSE, main = main, type = type, ylim = ylimcurv, ...)
        axis(side = 1, at = seq(0, 1, by = 0.1), labels = FALSE)
        axis(side = 2, labels = TRUE)
        # abline(h = 0, col = "red")
        plot(x = tortuosity_metrics$arclength/max(tortuosity_metrics$arclength, na.rm = T), tortuosity_metrics$torsion, ylab = ylabtors, axes = FALSE, type = type, ylim = ylimtors, ...)
        axis(side = 1, at = seq(0, 1, by = 0.1), labels = TRUE)
        axis(side = 2, labels = TRUE)
        mtext(text = "Normalized Arclength", line = -10, cex = 1.2)
        dev.off()
      }
      else if(format == "svg"){
        svg(filename = paste(filename, ".svg", sep = ""), pointsize = pointsize, width = width, height = height)
        par(mfrow = c(2, 1),     # 2x2 layout
            oma = c(5, 1, 1, 1), # 4 rows of text at the outer left and bottom margin
            mar = c(1, 3, 2, 0), # space for one row of text at ticks and to separate plots
            mgp = c(2, 1, 0),
            las = las) 
        plot(x = tortuosity_metrics$arclength/max(tortuosity_metrics$arclength, na.rm = T), tortuosity_metrics$curvature, xlab = NULL, ylab = ylabcurv, axes = FALSE, main = main, type = type, ylim = ylimcurv, ...)
        axis(side = 1, at = seq(0, 1, by = 0.1), labels = FALSE)
        axis(side = 2, labels = TRUE)
        abline(h = 0, col = "red")
        plot(x = tortuosity_metrics$arclength/max(tortuosity_metrics$arclength, na.rm = T), tortuosity_metrics$torsion, ylab = ylabtors, axes = FALSE, type = type, ylim = ylimtors, ...)
        axis(side = 1, at = seq(0, 1, by = 0.1), labels = TRUE)
        axis(side = 2, labels = TRUE)
        title(xlab = "Normalized Arclength", outer = TRUE)
        dev.off()
      }
    }else{
      par(mfrow = c(2, 1),     # 2x2 layout
          oma = c(5, 1, 1, 1), # 4 rows of text at the outer left and bottom margin
          mar = c(1, 3, 2, 0), # space for one row of text at ticks and to separate plots
          mgp = c(2, 1, 0),
          las = las) 
      plot(x = tortuosity_metrics$arclength/max(tortuosity_metrics$arclength, na.rm = T), tortuosity_metrics$curvature, xlab = NULL, ylab = ylabcurv, axes = FALSE, main = main, type = type, ylim = ylimcurv, las = las, ...)
      axis(side = 1, at = seq(0, 1, by = 0.1), labels = FALSE)
      axis(side = 2, labels = TRUE)
      abline(h = 0, col = "red")
      plot(x = tortuosity_metrics$arclength/max(tortuosity_metrics$arclength, na.rm = T), tortuosity_metrics$torsion, ylab = ylabtors, axes = FALSE, type = type, ylim = ylimtors, las = las, ...)
      axis(side = 1, at = seq(0, 1, by = 0.1), labels = TRUE)
      axis(side = 2, labels = TRUE)
      title(xlab = "Normalized Arclength", outer = TRUE)
    }
  }
  
  
  ## Plotting curvature and curvature check in one graph
  if(plot_type == 2){
    if(save_plot){
      png(filename = "curve_curvechk_vs_arclength.png", width = 8, height = 6, res = 300, units = "in")
      par(mfrow = c(2, 1),     # 2x2 layout
          oma = c(5, 1, 1, 1), # 4 rows of text at the outer left and bottom margin
          mar = c(1, 3, 1, 0), # space for one row of text at ticks and to separate plots
          mgp = c(2, 1, 0))
      plot(x = tortuosity_metrics$arclength/max(tortuosity_metrics$arclength, na.rm = T), tortuosity_metrics$curvature, xlab = NULL, ylab = "Curvature", axes = FALSE, main = "Curvature and Curvature Check vs Normalized Arclength", type = type)
      axis(side = 1, at = seq(0, 1, by = 0.1), labels = FALSE)
      axis(side = 2, labels = TRUE)
      plot(x = tortuosity_metrics$arclength/max(tortuosity_metrics$arclength, na.rm = T), tortuosity_metrics$curvature_check, ylab = "Curvature Check", axes = FALSE, type = type)
      axis(side = 1, at = seq(0, 1, by = 0.1), labels = TRUE)
      axis(side = 2, labels = TRUE)
      title(xlab = "Normalized Arclength", outer = TRUE)
      dev.off()
    }else{
      par(mfrow = c(2, 1),     # 2x2 layout
          oma = c(5, 1, 1, 1), # 4 rows of text at the outer left and bottom margin
          mar = c(1, 3, 1, 0), # space for one row of text at ticks and to separate plots
          mgp = c(2, 1, 0))
      plot(x = tortuosity_metrics$arclength/max(tortuosity_metrics$arclength, na.rm = T), tortuosity_metrics$curvature, xlab = NULL, ylab = "Curvature", axes = FALSE, main = "Curvature and Curvature Check vs Normalized Arclength", type = type)
      axis(side = 1, at = seq(0, 1, by = 0.1), labels = FALSE)
      axis(side = 2, labels = TRUE)
      plot(x = tortuosity_metrics$arclength/max(tortuosity_metrics$arclength, na.rm = T), tortuosity_metrics$curvature_check, ylab = "Curvature Check", axes = FALSE, type = type)
      axis(side = 1, at = seq(0, 1, by = 0.1), labels = TRUE)
      axis(side = 2, labels = TRUE)
      title(xlab = "Normalized Arclength", outer = TRUE)
    }
  }
  
  
  
  ## Plotting torsion and torsion check in one graph
  if(plot_type == 3){
    if(save_plot){
      png(filename = "torse_torsechk_vs_arclength.png", width = 8, height = 6, res = 300, units = "in")
      par(mfrow = c(2, 1),     # 2x2 layout
          oma = c(5, 1, 1, 1), # 4 rows of text at the outer left and bottom margin
          mar = c(1, 3, 1, 0), # space for one row of text at ticks and to separate plots
          mgp = c(2, 1, 0))
      plot(x = tortuosity_metrics$arclength/max(tortuosity_metrics$arclength, na.rm = T), tortuosity_metrics$torsion, xlab = NULL, ylab = "Torsion", axes = FALSE, main = "Torsion and Torsion Check vs Normalized Arclength", type = type)
      axis(side = 1, at = seq(0, 1, by = 0.1), labels = FALSE)
      axis(side = 2, labels = TRUE)
      plot(x = tortuosity_metrics$arclength/max(tortuosity_metrics$arclength, na.rm = T), tortuosity_metrics$torsion_check, ylab = "Torsion Check", axes = FALSE, type = type)
      axis(side = 1, at = seq(0, 1, by = 0.1), labels = TRUE)
      axis(side = 2, labels = TRUE)
      title(xlab = "Normalized Arclength", outer = TRUE)
      dev.off()
    }else{
      par(mfrow = c(2, 1),     # 2x2 layout
          oma = c(5, 1, 1, 1), # 4 rows of text at the outer left and bottom margin
          mar = c(1, 3, 1, 0), # space for one row of text at ticks and to separate plots
          mgp = c(2, 1, 0))
      plot(x = tortuosity_metrics$arclength/max(tortuosity_metrics$arclength, na.rm = T), tortuosity_metrics$torsion, xlab = NULL, ylab = "Torsion", axes = FALSE, main = "Torsion and Torsion Check vs Normalized Arclength", type = type)
      axis(side = 1, at = seq(0, 1, by = 0.1), labels = FALSE)
      axis(side = 2, labels = TRUE)
      plot(x = tortuosity_metrics$arclength/max(tortuosity_metrics$arclength, na.rm = T), tortuosity_metrics$torsion_check, ylab = "Torsion Check", axes = FALSE, type = type)
      axis(side = 1, at = seq(0, 1, by = 0.1), labels = TRUE)
      axis(side = 2, labels = TRUE)
      title(xlab = "Normalized Arclength", outer = TRUE)
    }
  }
}