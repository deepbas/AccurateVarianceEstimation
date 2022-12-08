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
p <- 5
rho <- seq(0.80, 0.95, by = 0.01)
n <- 2e4
nrep <- 100
omega <- diag(p)
#%-------------------------------------------------

true_Sigmas  <- list(length = length(rho))
Vs  		 <- list(length = length(rho))
phis		 <- list(length = length(rho))
sims_for_rho <- list(length = length(rho))
# generating VAR(1) process

detectCores()

# leave 2 cores free for computer to work well
# this way you can walk Netflix and run code!

# this tells the machine to register the 4 cores
registerDoParallel(cores = detectCores()-2)

# for all values of rho
for(s in 1:length(rho))
{
	print(paste("Doing rho = ", rho[s]))
	temp 				<- sigphi(p, rho[s])
	phis[[s]] 			<- temp[[1]]
	true_Sigmas[[s]] 	<- temp[[2]]
	Vs[[s]] 			<- temp[[3]]


	exact_autcov <- true_autocov(phi = phis[[s]], V = Vs[[s]], lag = n-1)
	## a doParallel for reps
	sims_for_rho[[s]] 	<- foreach(st = 1:nrep) %dopar% 
	{
		chain <- as.matrix(mAr.sim(rep(0,p), as.matrix(phis[[s]]), omega, N = n))
		est_var(chain = as.matrix(chain[,1]), phi = phis[[s]], Sigma = true_Sigmas[[s]], exact_autcov = exact_autcov)
	}	
}

save(file = "var_batch", sims_for_rho, phis, true_Sigmas, rho)


