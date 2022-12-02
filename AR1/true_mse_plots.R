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
mse_plot(n, rho = .99, b.seq = b.seq)




