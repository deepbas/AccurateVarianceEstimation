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
rho <- seq(0.85, 0.99, by = 0.01)
n <- 2e4
nrep <- 2
omega <- diag(p)
#%-------------------------------------------------

true_Sigmas  <- list(length = length(rho))
phis		 <- list(length = length(rho))
sims_for_rho <- list(length = length(rho))
# generating VAR(1) process

detectCores()

# leave 2 cores free for computer to work well
# this way you can walk Netflix and run code!

# this tells the machine to register the 4 cores
registerDoParallel(cores = detectCores()-2)

ar1 <- function(N, phi, omega, start)
{
	out <- numeric(length = N)
	out[1] <- start
	eps <- rnorm(N, sd = sqrt(omega))
	for(t in 2:N)
	{
		out[t] <- phi*out[t-1] + eps[t]
	}
	return(out)
}
# for all values of rho
for(s in 1:length(rho))
{
	print(paste("Doing rho = ", rho[s]))
	temp 	<- sigphi(p, rho[s])
	phis[[s]] 	<- temp[[1]]
	true_Sigmas[[s]] 	<- temp[[2]]

	## a doParallel for reps
	sims_for_rho[[s]] 	<- foreach(st = 1:nrep) %dopar% 
	{

		chain <- as.matrix(ar1(N = n, phi = phis[[s]][1,1], omega = omega[1,1], start = 0))
		est_var(chain = chain, phi = phis[[s]], Sigma = true_Sigmas[[s]])
	}	
}

save(file = "ar1_out", sims_for_rho, phis, true_Sigmas, rho)


