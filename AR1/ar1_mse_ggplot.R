set.seed(16235, kind = "L'Ecuyer-CMRG" )
source("batchSizes.R")
source("../backFuncs.R")
#%-------------------------------------------------
pkgs <- c("doParallel", "Matrix", "ts.extend", "mAr", "mcmcse")

 if(sum(as.numeric(!pkgs %in% installed.packages())) != 0) {
    installer <- pkgs[!pkgs %in% installed.packages()]
    for(i in 1:length(installer)) {
      install.packages(installer, dependencies = T)
      break()}
    sapply(pkgs, require, character = T)
  } else {
    sapply(pkgs, require, character = T)
  }

#%-------------------------------------------------
# Simulation settings
p <- 1
rho <- .99
n <- 1e3
nrep <- 1 # increase for estimation
omega <- diag(p)
#%-------------------------------------------------

mse_plot <- function(n = 1e3, rho = .99, b.seq)
{
	Sigma <- 1/(1 - rho)^2
	exact.autocov <- rho^(0:(n-1))*(1/(1 - rho^2))
	gamma.exact <- -2*sum( (0:(n-1)* exact.autocov))

	true.mse <- matrix(0, ncol = length(b.seq), nrow = 3)
	c <- 1/2
	for(r in 1:3)
	{
		for(i in 1:length(b.seq))
		{
			true.mse[r, i] <- funBMexact(b = b.seq[i], x = exact.autocov, 
				y = Sigma, r = r, c = c)
		}
	}


	b.bm.exact <- sapply(1:3, function(k) optim(par = c(40),
	                     fn = funBMexact, 
	                     x  = exact.autocov, 
	                     y  = Sigma, 
	                     r  = k, c = c, 
	                     method = "Brent", lower = 5, 
	                     upper = n/2)$par)	

	# Estimated batch size - our method
	b.bm <- sapply(1:3, function(k) optim(par = c(40),
	               fn = funBMi,  
	               x  = exact.autocov, # do ar.autocov for estimated quantities
	               y  = Sigma, # do Sigma.pilot for estimated quantities
	               r  = k, c = c, 
	               method = "Brent", lower = 5, 
	               upper = n/2)$par)	

	# Current first order method
	b.curr.bm <- sapply(1:3, function(k)  optim(par = c(40), 
	                fn = funCurrbm, 
	                n  = n, 
	                x  = gamma.exact, # gamma.pilot for estim
	                y  = Sigma, # Sigma.pilot for estimated quantities
	                r  = k, c = c, 
	                method = "Brent", lower = 5, 
	                upper = n/2 )$par)	
	for(r in 1:3)
	{
	pdf(paste("plots/mse_r", r, "n",n,".pdf", sep = ""), height = 5, width = 5)
	  
	plot(b.seq, true.mse[r, ], type = 'l', 
		ylab = "Meas Squared Error", xlab = "Batch Size")
	points(b.bm.exact[r], funBMexact(b = b.bm.exact[r], x = exact.autocov, 
				y = Sigma, r = r, c = c), col = 1, pch = 7, lwd = 2)
	points(b.bm[r], funBMexact(b = b.bm[r], x = exact.autocov, 
				y = Sigma, r = r, c = c), col = 2, pch = 8, lwd = 2)
	points(b.curr.bm[r], funBMexact(b = b.curr.bm[r], x = exact.autocov, 
				y = Sigma, r = r, c = c), col = 3, pch = 13, lwd = 2)
	legend("top", legend = c("Exact", "Higher-order", "First-Order"),
		pch = c(7,8,13), col = 1:3, bty = "n")
	dev.off()
	}
}

n <- 1e3
b.seq <- seq(5, n/4, by = 5)
mse_plot(n, rho = .99, b.seq = b.seq)

n <- 1e4
b.seq <- seq(5, n/10, by = 5)



# Function to handle shared legend across multiple plots
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = "bottom"))
  legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  
  main_plot <- do.call(arrangeGrob, c(plots, ncol = length(plots)))
  
  combined <- arrangeGrob(main_plot, legend, 
                          heights = unit.c(unit(1, "npc") - lheight, lheight),
                          ncol = 1)
  return(combined)
}


mse_plot_new <- function(n = 1e3, rho = .99, b.seq)
{
  Sigma <- 1/(1 - rho)^2
  exact.autocov <- rho^(0:(n-1))*(1/(1 - rho^2))
  gamma.exact <- -2*sum( (0:(n-1)* exact.autocov))
  
  true.mse <- matrix(0, ncol = length(b.seq), nrow = 3)
  c <- 1/2
  for(r in 1:3)
  {
    for(i in 1:length(b.seq))
    {
      true.mse[r, i] <- funBMexact(b = b.seq[i], x = exact.autocov, 
                                   y = Sigma, r = r, c = c)
    }
  }
  
  
  b.bm.exact <- sapply(1:3, function(k) optim(par = c(40),
                                              fn = funBMexact, 
                                              x  = exact.autocov, 
                                              y  = Sigma, 
                                              r  = k, c = c, 
                                              method = "Brent", lower = 5, 
                                              upper = n/2)$par)	
  
  # Estimated batch size - our method
  b.bm <- sapply(1:3, function(k) optim(par = c(40),
                                        fn = funBMi,  
                                        x  = exact.autocov, # do ar.autocov for estimated quantities
                                        y  = Sigma, # do Sigma.pilot for estimated quantities
                                        r  = k, c = c, 
                                        method = "Brent", lower = 5, 
                                        upper = n/2)$par)	
  
  # Current first order method
  b.curr.bm <- sapply(1:3, function(k)  optim(par = c(40), 
                                              fn = funCurrbm, 
                                              n  = n, 
                                              x  = gamma.exact, # gamma.pilot for estim
                                              y  = Sigma, # Sigma.pilot for estimated quantities
                                              r  = k, c = c, 
                                              method = "Brent", lower = 5, 
                                              upper = n/2 )$par)	
  
  # Data frame for MSE curves
  
 
  
  plot_list <- list()
  
  for(r in 1:3) {
    
    # Data frame for MSE curves
    df_mse <- data.frame(BatchSize = b.seq,
                         MSE = true.mse[r, ])
    
    # Data frame for the points
    df_points <- data.frame(BatchSize = c(b.bm.exact[r], b.bm[r], b.curr.bm[r]),
                            MSE = sapply(c(b.bm.exact[r], b.bm[r], b.curr.bm[r]), function(b) {
                              funBMexact(b = b, x = exact.autocov, y = Sigma, r = r, c = c)
                            }),
                            Method = factor(c("Exact", "Higher-order", "First-Order")),
                            Shape = c(16, 17, 18))
    
    # Create ggplot
    p <- ggplot(df_mse, aes(x = BatchSize, y = MSE)) +
      geom_line() +
      geom_point(data = df_points, aes(x = BatchSize, y = MSE, color = Method, shape = Method), size = 4) +
      scale_color_manual(values = c("Exact" = 1, "Higher-order" = 2, "First-Order" = 3)) +
      scale_shape_manual(values = c("Exact" = 16, "Higher-order" = 17, "First-Order" = 18)) +
      labs(
        x = "b",
        y = "MSE",
        subtitle = paste("r =", r)
      ) + #ylim(c(0.1e7, 1e8)) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(angle = 45, size = 8),
        legend.position = "none",
        plot.subtitle = element_text(size = 10, hjust = 0.5)
      )
    plot_list[[r]] <- p
  }
  
  # Generate the combined grob
  combined_plot <- grid_arrange_shared_legend(plot_list[[1]], plot_list[[2]], plot_list[[3]])
  
  # Convert to a ggplot object using cowplot's ggdraw and as.grob
  combined_gg <- ggdraw() + draw_grob(combined_plot)
  
  return(combined_gg)
}



n <- 1e3
b.seq <- seq(5, n/4, by = 5)
plot_list <- mse_plot_new(n, rho = .95, b.seq = b.seq)


# Use ggsave to save the ggplot object
ggsave("plots/MSEAR1_1e3rho95.pdf", plot = plot_list, width = 11, height = 5)
