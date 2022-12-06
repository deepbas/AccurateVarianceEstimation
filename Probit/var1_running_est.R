set.seed(16235, kind = "L'Ecuyer-CMRG" )
source("batchSizes.R")
source("../backFuncs.R")


#%-------------------------------------------------
pkgs <- c("doParallel", "Matrix", "ts.extend", "mAr", "mcmcse", "TruncatedNormal", "mvtnorm")

 if(sum(as.numeric(!pkgs %in% installed.packages())) != 0) {
    installer <- pkgs[!pkgs %in% installed.packages()]
    for(i in 1:length(installer)) {
      install.packages(installer, dependencies = T)
      break()}
    sapply(pkgs, require, character = T)
  } else {
    sapply(pkgs, require, character = T)
  }


## -----------------------------------------------------
# Function to simulate the posterior samples in the lupus example

b.gen <- function(n.sim){
  data("lupus")
  Y <- lupus[,1]; # response data
  X <- as.matrix(lupus[,-1])  # construct design matrix
  n <- nrow(X)
  p <- ncol(X)
  
  #X <- diag(2*Y-1) %*% X; # incorporate response into design matrix
  b <- matrix(NA, n.sim, 3)
  
  # initialize
  b[1,] <- c(-1.778, 4.374, 2.482)
  for (j in 2:n.sim){
    m <- X%*%b[j-1,]
    lower <- ifelse(Y==1, 0, -Inf)
    upper <- ifelse(Y==1, Inf, 0)
    z <- sapply(1:n, function(i) rtnorm(1,mu= m[i], sd = 1, lb = lower[i], ub = upper[i]))
    shape <- 55/2
    foo <- sum((z - X%*%solve(t(X)%*%X)%*%t(X)%*%z)^2)
    scale <- 1/2 * foo
    
    g.sq <- rgamma(1, shape, scale)
    
    z.prime <- sqrt(g.sq)*z
    
    b.prime <- rmvnorm(n = 1, mean = solve(t(X)%*%X)%*%t(X)%*%z.prime, sigma = solve(t(X)%*%X))
    b[j,] <- b.prime
  }
  return(b)
}




#%-------------------------------------------------
############################################
## file runs running estimation
## of optimal batch sizes and resulting variance
## estimators
############################################


# nseq = sequence of ns where to calculate batch size etc
running_est <- function(chain, nseq)
{
	
	est_n <- list(length(nseq))

	for(nind in 1:length(nseq))
	{
		print(nind)
		n.ind <- nseq[nind]
		sub.chain <- as.matrix(chain[1:n.ind, ])
		batchObj <- batch_sizes(sub.chain)

		secondBM 	<- lapply(1:3, function(k) mcse.multi(sub.chain, size = ceiling(batchObj[[1]][k,1]), r = k)$cov)
		firstBM 	<- lapply(1:3, function(k) mcse.multi(sub.chain, size = ceiling(batchObj[[1]][k,2]), r = k)$cov)
		politBM 	<- mcse.multi(chain, size = ceiling(batchObj[[2]]), r = 2)$cov

		est_n[[nind]] <- list("batches" = batchObj,
		"secondBM" = secondBM,
		"firstBM" = firstBM, "politBM" = politBM)
	}
	return(est_n)
}


# Simulation settings
nsim <- 1e5
nrep <- 2

# generating VAR(1) process

detectCores()

# leave 2 cores free for computer to work well
# this way you can walk Netflix and run code!

# this tells the machine to register the 4 cores
registerDoParallel(cores = detectCores()-2)


# calculating phi and truths

# where to find batch sizes etc
nseq <- floor(seq(1e3, 1e5, length = 250))


## a doParallel for reps

sims_for_lupus 	<- foreach(st = 1:nrep) %dopar% 
{
	print(st)
	lupus_chain <- as.matrix(b.gen(nsim))
	running_est(chain = lupus_chain, nseq = nseq)
	
}	


save(file = "var1_running", sims_for_lupus, nseq)
