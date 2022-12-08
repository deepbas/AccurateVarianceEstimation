set.seed(16235, kind = "L'Ecuyer-CMRG" )
source("batchSizes.R")
source("probit_chain.R")
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

		secondBM 	<- lapply(1:3, function(k) mcse.multi(sub.chain, size = ceiling(batchObj[k,1]), r = k)$cov)
		firstBM 	<- lapply(1:3, function(k) mcse.multi(sub.chain, size = ceiling(batchObj[k,2]), r = k)$cov)

		est_n[[nind]] <- list("batches" = batchObj,
		"secondBM" = secondBM,
		"firstBM" = firstBM)
	}
	return(est_n)
}


# Simulation settings
nsim <- 1e5
nrep <- 100

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


save(file = "probit_running", sims_for_lupus, nseq)
