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
############################################


# n = length of Markov chain
# p = dimension of the problem
# Sigma =  true matrix to be estimated
# phi = phi in VAR(1)

running_est <- function(chain, phi, Sigma, 
	nseq = c(1e3, 3e3, 5e3, 1e4, 2e4, 5e4))
{
	# Find all batch sizes
	est_n <- list(length(nseq))
	for(nind in 1:length(nseq))
	{
		n.ind <- nseq[nind]
		print(n.ind)

		sub.chain <- matrix(chain[1:n.ind, ])
		batchObj <- batch_sizes(sub.chain, phi, Sigma)

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
p <- 1
rho <- .99
nrep <- 5
omega <- diag(p)
#%-------------------------------------------------

true_Sigmas  <- list(length = length(rho))
phis		 <- list(length = length(rho))
sims_for_n <- list(length = length(rho))
# generating VAR(1) process

detectCores()

# leave 2 cores free for computer to work well
# this way you can walk Netflix and run code!

# this tells the machine to register the 4 cores
registerDoParallel(cores = detectCores()-2)


# for all values of rho
phis 	<- rho
true_Sigmas	<- 1/(1 - rho)^2
nseq <- floor(seq(1e3, 5e4, length = 300))

## a doParallel for reps
sims_for_n 	<- foreach(st = 1:nrep) %dopar% 
{
	print(st)
	chain <- as.matrix(ar1(N = max(nseq), phi = phis, omega = omega[1,1], start = 0))
	running_est(chain = chain, phi = as.matrix(phis), Sigma = as.matrix(true_Sigmas), nseq = nseq)
}	


save(file = "ar1_running", sims_for_n, phis, true_Sigmas, nseq)
