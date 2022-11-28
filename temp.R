library(mcmcse)
n <- 1e4
b <- 50
rho <- .99

reps <- 1e3
truth <- 1/(1 - rho)^2

est <- matrix(0, nrow = reps, ncol = 2)
for(r in 1:reps)
{
  out <- numeric(length = n)
  out[1] <- 0
  err <- rnorm(n)
  for(i in 2:n)
  {
    out[i] <- rho*out[i-1] + err[i]
  }

  est[r, 1] <- as.numeric(mcse.multi(out, r = 1,size = b)$cov)
  est[r, 2] <- as.numeric(mcse.multi(out, r = 2,size = b)$cov)
}

### checking covariances
cov(est)[1, ]
cov.r <- 1:5
(2*truth^2*b/n)*(1/cov.r + b/cov.r* 1 / (n*cov.r + b))

#### checking variance estimates
foo <- apply(est, 2, var)



r <- 1
c <- 1/2
# theory <-  (2*truth^2*b/n)*( (b*c^2/(1-c)^2)/(n*r-b) - 
#     (2*b*c)/((1-c)^2*r*(n*r-b)) + b/((1-c)^2*(n-b)) + 
#     c^2/(1-c)^2 - 2*c/((1-c)^2*r) + 1/(1-c)^2 )

theory <- (2*truth^2*b/n) * (1 + b/(n-b) + 
  c*(c-2)*( (1/r) + b/(r*(n*r - b))) )/(1 - c)^2
theory/var(est[,1])

r <- 2
c <- 1/2
# theory <-  (2*truth^2*b/n)*( (b*c^2/(1-c)^2)/(n*r-b) - 
#     (2*b*c)/((1-c)^2*r*(n*r-b)) + b/((1-c)^2*(n-b)) + 
#     c^2/(1-c)^2 - 2*c/((1-c)^2*r) + 1/(1-c)^2 )

theory <- (2*truth^2*b/n)* (1 + b/(n-b) + 
  c*(c-2)*( (1/r) + b/(r*(n*r - b))) ) /(1 - c)^2 


theory/var(est[,2])


#################
## Checking bias
n <- 1e4
b <- 100
c <- 1/2
r <- 1
rho <- .993

temp <- function(b)
{
  r <- 1
  br <- b/r
  a <- n/b
  ar <- n/br
  phi <- rho
  x <- phi^(0:(n-1))*(1/(1 - phi^2))
  y <- 1/(1-rho)^2

  gamma0b1 <- 2*sum(x[2:b])
  gamma0b1r <- 2*sum(x[2:br])
  gamma0n1 <- 2*sum(x[(b+1):n])
  gamma0n1r <- 2*sum(x[(br+1):n])

  gamma1b1 <- 2*sum((1:(b-1))*x[2:b])
  gamma1b1r <- 2*sum((1:(br-1))*x[2:br])
  gamma1n1 <- 2*sum(x[(b+1):n]*(b:(n-1)))
  gamma1n1r <- 2*sum(x[(br+1):n]*(br:(n-1)))

  bias1 <- x[1] + gamma0b1 - (a + 1)/(a*b)*gamma1b1 - 
          1/(a - 1)*(gamma0n1  -  gamma1n1/n )
  bias2 <- x[1] + gamma0b1r - (ar + 1)/(ar*br)*gamma1b1r - 
          1/(ar - 1)*(gamma0n1r  -  gamma1n1r/n )          
  theory <- y - ( (1/(1-c))*bias1 - (c/(1-c))*bias2)
  theory <- -theory

}

temp2 <- function(b)
{
  a <- n/b
    phi <- rho
  x <- phi^(0:(n-1))*(1/(1 - phi^2))
  x <- -2*sum( (0:(n-1)* x))
  out <- (x*(a+1)/n)*(1-r*c)/(1-c)
}

bs <- seq(10, n/10, length = 1e3)
foo <- numeric(length = length(bs))
foo2 <- numeric(length = length(bs))
for(i in 1:length(bs))
{
  foo[i] <- temp(bs[i])
  foo2[i] <- temp2(bs[i])
}
plot(bs, foo, type = "l", ylim = range(c(foo, foo2)))
lines(bs, foo2, col = "red")
# theory - (truth - mean(est[,1]))



temp <- function(b)
{
  (2*y^2*b/n)* (1 + b/(n-b) + 
  c*(c-2)*( (1/r) + b/(r*(n*r - b))) ) /(1 - c)^2 
}
temp2 <- function(b)
{
  2*b/n*(y^2)*(1/r + (r-1)/(r*(1-c)^2)) 
}
bs <- seq(10, n/2, length = 1e3)
foo <- numeric(length = length(bs))
foo2 <- numeric(length = length(bs))
for(i in 1:length(bs))
{
  foo[i] <- temp(bs[i])
  foo2[i] <- temp2(bs[i])
}
plot(bs, foo, type = "l", ylim = range(c(foo, foo2)))
lines(bs, foo2, col = "red")





temp <- function(b)
{
  br <- b/r
  a <- n/b
  ar <- n/br
  phi <- rho
  x <- phi^(0:(n-1))*(1/(1 - phi^2))
  y <- 1/(1-rho)^2

  gamma0b1 <- 2*sum(x[2:b])
  gamma0b1r <- 2*sum(x[2:br])
  gamma0n1 <- 2*sum(x[(b+1):n])
  gamma0n1r <- 2*sum(x[(br+1):n])

  gamma1b1 <- 2*sum((1:(b-1))*x[2:b])
  gamma1b1r <- 2*sum((1:(br-1))*x[2:br])
  gamma1n1 <- 2*sum(x[(b+1):n]*(b:(n-1)))
  gamma1n1r <- 2*sum(x[(br+1):n]*(br:(n-1)))

  bias1 <- x[1] + gamma0b1 - (a + 1)/(a*b)*gamma1b1 - 
          1/(a - 1)*(gamma0n1  -  gamma1n1/n )
  bias2 <- x[1] + gamma0b1r - (ar + 1)/(ar*br)*gamma1b1r - 
          1/(ar - 1)*(gamma0n1r  -  gamma1n1r/n )          
  theory <- y - ( (1/(1-c))*bias1 - (c/(1-c))*bias2)
  bias <- -theory

  variance <- (2*y^2*b/n)* (1 + b/(n-b) + 
  c*(c-2)*( (1/r) + b/(r*(n*r - b))) ) /(1 - c)^2 

  bias^2 + variance
}

temp2 <- function(b)
{
    a <- n/b
    phi <- rho
  x <- phi^(0:(n-1))*(1/(1 - phi^2))
  x <- -2*sum( (0:(n-1)* x))
  y <- 1/(1-rho)^2
  bias <- (x*(a+1)/n)*(1-r*c)/(1-c)
  variance <- 2*b/n*(y^2)*(1/r + (r-1)/(r*(1-c)^2)) 
  bias^2 + variance
}

bs <- seq(5, n/5, length = 5e2)
foo <- numeric(length = length(bs))
foo2 <- numeric(length = length(bs))
for(i in 1:length(bs))
{
  foo[i] <- temp(bs[i])
  foo2[i] <- temp2(bs[i])
}
plot(bs, foo, type = "l", ylim = range(c(foo, foo2)))
lines(bs, foo2, col = "red")

