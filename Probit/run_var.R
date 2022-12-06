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
