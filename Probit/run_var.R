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

#%-------------------------------------------------
# Simulation settings
nsim <- 2e4
nrep <- 2

detectCores()

# leave 2 cores free for computer to work well
# this way you can walk Netflix and run code!

# this tells the machine to register the 4 cores
registerDoParallel(cores = detectCores()-2)


sims_for_lupus <- foreach(st = 1:nrep) %dopar% {
lupus_chain <- as.matrix(b.gen(nsim))
est_var(chain = lupus_chain)
}

save(file = "var_out", sims_for_lupus)
