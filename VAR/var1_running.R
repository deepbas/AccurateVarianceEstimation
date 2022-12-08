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
############################################
## file runs running estimation
## of optimal batch sizes and resulting variance
## estimators
############################################


# n = length of Markov chain
# p = dimension of the problem
# Sigma =  true matrix to be estimated
# phi = phi in VAR(1)
# nseq = sequence of ns where to calculate batch size etc
running_est <- function(chain, phi, Sigma, nseq)
{
	
	est_n <- list(length(nseq))

	for(nind in 1:length(nseq))
	{
		print(nind)
		n.ind <- nseq[nind]
		sub.chain <- as.matrix(chain[1:n.ind, ])
		batchObj <- batch_sizes(sub.chain, phi, Sigma, exact_autcov[1:n.ind, ])

		exactBM		<- lapply(1:3, function(k) mcse.multi(sub.chain, size = ceiling(batchObj[k,1]), r = k)$cov)
		secondBM 	<- lapply(1:3, function(k) mcse.multi(sub.chain, size = ceiling(batchObj[k,2]), r = k)$cov)
		firstBM 	<- lapply(1:3, function(k) mcse.multi(sub.chain, size = ceiling(batchObj[k,3]), r = k)$cov)
		
		est_n[[nind]] <- list("batches" = batchObj,
		"exactBM" = exactBM, "secondBM" = secondBM,
		"firstBM" = firstBM)
	}
	return(est_n)
}


# Simulation settings
p <- 5
rho <- .95
nrep <- 100
omega <- diag(p)
#%-------------------------------------------------

true_Sigmas  <- list(length = length(rho))
Vs  		 <- list(length = length(rho))
phis		 <- list(length = length(rho))
sims_for_n 	 <- list(length = length(rho))
#%-------------------------------------------------

# generating VAR(1) process

detectCores()

# leave 2 cores free for computer to work well
# this way you can walk Netflix and run code!

# this tells the machine to register the 4 cores
registerDoParallel(cores = detectCores()-2)


# calculating phi and truths
temp 				<- sigphi(p, rho)
phis				<- temp[[1]]
true_Sigmas <- temp[[2]]
Vs					<- temp[[3]]

# where to find batch sizes etc
nseq <- floor(seq(1e3, 1e5, length = 250))

# exact autocov for each function
exact_autcov <- true_autocov(phi = phis, V = Vs, lag = max(nseq)-1)

## a doParallel for reps
sims_for_n 	<- foreach(st = 1:nrep) %dopar% 
{
	print(st)
	chain <- as.matrix(mAr.sim(rep(0,p), as.matrix(phis), omega, N = max(nseq)))
	running_est(chain = chain, phi = phis, Sigma = true_Sigmas, nseq = nseq)
}	


save(file = "var1_running", sims_for_n, phis, true_Sigmas, nseq)
