############################################
## file runs one loop of the batching
############################################


# n = length of Markov chain
# p = dimension of the problem
# Sigma =  true matrix to be estimated
# phi = phi in VAR(1)
batch_sizes <- function(chain, phi, Sigma)
{
	p 			<- dim(chain)[2]
	n 			<- dim(chain)[1]
	c  			<-  1/2
	b 			<- seq(1,n^(0.62),1)
	a 			<- n/b
 	ubound 		<- 2*sqrt(log(n)/n) # needed for lag-based

	m 			<- numeric(length = p)
	phi.i 		<- list(length = p)
	sigma.e 	<- numeric(length = p)
	Sigma.pilot <- numeric(length = p)
	gamma.pilot <- numeric(length = p)
	
	b.politis	<- numeric(length = p)
	b.bm 		<- matrix(0, nrow = p, ncol = 3, 
					dimnames = list(NULL, c("r = 1", "r = 2", "r = 3")))
	b.bm.exact	<- matrix(0, nrow = p, ncol = 3, 
				dimnames = list(NULL, c("r = 1", "r = 2", "r = 3")))
	b.curr.bm	<- matrix(0, nrow = p, ncol = 3, 
				dimnames = list(NULL, c("r = 1", "r = 2", "r = 3")))
	
	# Loop for each component separately
	for(i in 1:p)
	{
		## chunk for calculating pilot gamma
		## some other approximations are also made on the 
		## way that will be useful later on

		#%--------------------------------
		ar.chain 	 <- ar(chain[,i], order.max = NULL, method = "yw")
		m		     <- ar.chain$order
		phi.i	  	 <- ar.chain$ar
		sigma.e	     <- ar.chain$var.pred
		Sigma.pilot  <- sigma.e/(1 - sum(phi.i))^2

		ar.autocovar <- acf(chain[,i], lag.max = n-1, 
						type = "covariance", plot = FALSE)$acf

		foo <- 0
		for(j in 1:m){
			for(k in 1:j){
			  foo <- foo + phi.i[j]*k*ar.autocovar[abs(k-j)+1]
			}
		}

  		gamma.pilot <- -2*(foo + (Sigma.pilot - ar.autocovar[1])/2 *
  			sum(1:m * phi.i)  )/(1 - sum(phi.i))
  		#%--------------------------------	


  		# Calculating component-wise batch size
  		# for each lugsail parameter
  		# will combine outside loop
		#%--------------------------------	

  		exact.autocov <- phi[1,1]^(0:(n-1))*(1/(1 - phi[1,1]^2))
		ar.autocovar   <- ARMA.autocov(n = n, ar = phi.i, 
							ma = 0, corr = FALSE)

		# dum <- -2*sum( (0:(n-1)* exact.autocov))
		# Exact theoretically optimal batch size

		# foo <- optim(par = c(40),
  #                    fn = funBMexact, 
  #                    x  = exact.autocov, 
  #                    y  = diag(Sigma)[i], 
  #                    r  = 2, c = c, 
  #                    method = "Brent",
  #                    lower = c(0), upper=c(5000))$par

		b.bm.exact[i, ] <- sapply(1:3, function(k) optim(par = c(100),
		                     fn = funBMexact, 
		                     x  = exact.autocov, 
		                     y  = diag(Sigma)[i], 
		                     r  = k, c = c, 
		                     method = "Brent", lower = 5, 
		                     upper = n/2)$par)	

		# Estimated batch size - our method
		b.bm[i, ] <- sapply(1:3, function(k) optim(par = c(40),
                       fn =funBMi,  
                       x  = ar.autocovar , 
                       y  = Sigma.pilot, 
                       r  = k, c = c, 
                       method = "Brent", lower = 5, 
		               upper = n/2)$par)	

		# Current first order method
		b.curr.bm[i, ] <- sapply(1:3, function(k)  optim(par = c(40), 
                        fn = funCurrbm,  
                        x  = gamma.pilot, 
                        y  = Sigma.pilot, 
                        r  = k, c = c, 
                        method = "Brent", lower = 5, 
		                upper = n/2 )$par)	

		# Politis method
		ar.autocorr <- abs(acf(chain[,i], lag.max = n-1, 
						type = "correlation", plot = FALSE)$acf)

		politis <- TRUE
		ind <- 0
		while(politis)
		{
			foo_again <- which(ar.autocorr < ubound)
			pol.r <-  foo_again[which(foo_again > ind)] [1]
			if( all(ar.autocorr[pol.r: (pol.r + 5)] < ubound) )
			{
				b.politis[i] <- 2*pol.r
				politis <- FALSE
			} else{
				ind <- pol.r + 5
			}
		}
		#%--------------------------------		
	}

	# Taking geometric means
	b.exact <- apply(b.bm.exact, 2, function(x) exp(mean(log(x))))
	b.bm.opt <- apply(b.bm, 2, function(x) exp(mean(log(x))))
	b.currbm.opt <- apply(b.curr.bm, 2, function(x) exp(mean(log(x))))
	b.politis <- exp(mean(log(b.politis)))

	b.sizes <- cbind(b.exact, b.bm.opt, b.currbm.opt)
	colnames(b.sizes) <- c("exact", "secondOpt", "firstOpt")

	return(list(b.sizes, b.politis))
}

est_var <- function(chain, phi, Sigma)
{
	# Find all batch sizes
	batchObj <- batch_sizes(chain, phi, Sigma)

	exactBM		<- lapply(1:3, function(k) mcse.multi(chain, size = ceiling(batchObj[[1]][k,1]), r = k)$cov)
	secondBM 	<- lapply(1:3, function(k) mcse.multi(chain, size = ceiling(batchObj[[1]][k,2]), r = k)$cov)
	firstBM 	<- lapply(1:3, function(k) mcse.multi(chain, size = ceiling(batchObj[[1]][k,3]), r = k)$cov)
	politBM 	<- mcse.multi(chain, size = ceiling(batchObj[[2]]), r = 2)$cov

	return(list("batches" = batchObj,
		"exactBM" = exactBM, "secondBM" = secondBM,
		"firstBM" = firstBM, "politBM" = politBM))
}




