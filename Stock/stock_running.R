set.seed(16235, kind = "L'Ecuyer-CMRG" )
source("batchSizes.R")
source("stock_chain.R")
source("../backFuncs.R")


#%-------------------------------------------------
pkgs <- c("doParallel", "Matrix", "ts.extend", "mAr", "mcmcse", "TruncatedNormal", "mvtnorm", "yfR")

 if(sum(as.numeric(!pkgs %in% installed.packages())) != 0) {
    installer <- pkgs[!pkgs %in% installed.packages()]
    for(i in 1:length(installer)) {
      install.packages(installer, dependencies = T)
      break()}
    sapply(pkgs, require, character = T)
  } else {
    sapply(pkgs, require, character = T)
  }



my_tickers <- c("AAPL", "MSFT", "TXN", "AMD")
first_date <- "2010-01-01"
last_date <- "2023-01-01"

# fetch data
df_yf <- yf_get(tickers = my_tickers, 
                first_date = first_date,
                last_date = last_date,
                freq_data = 'daily')

df_wide <- yf_convert_to_wide(df_yf)

my_data <- df_wide$price_close[, -1]
y <- scale(my_data, center = FALSE, scale = apply(my_data, 2, sd))
R <- cov2cor(cov(y))
V <- cov(my_data)
D <- diag(diag(V))

tau2 <- 1
alpha <- 1
beta <- 0.90
v2 <- .005
nsim <- 1e5
nrep <- 1
#h <- c(.0023, 0.001, 0.001, 0.0024, 0.002, .0015)
h <- c(.0025, 0.0020, 0.0018, 0.0024, 0.002, .0025)


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
nrep <- 2

# generating VAR(1) process

detectCores()

# leave 2 cores free for computer to work well
# this way you can walk Netflix and run code!

# this tells the machine to register the 4 cores
registerDoParallel(cores = detectCores()-2)


# calculating phi and truths

# where to find batch sizes etc
nseq <- floor(seq(6e4, 1e5, length = 5))

## a doParallel for reps

sims_for_stocks 	<- foreach(st = 1:nrep) %dopar% 
{
	print(st)
  
  chain_crude <- CorMCMC(y = y, tau2 = tau2, alpha = alpha, beta = beta, port = .90, nsim = nsim, h = h)
  stock_chain <- as.matrix(do.call(rbind, lapply(1:nsim, function(idx) lower(chain_crude[[1]][idx, 1:6], D))))
	running_est(chain = stock_chain, nseq = nseq)
}	


save(file = "stock_running", sims_for_stocks, nseq)
