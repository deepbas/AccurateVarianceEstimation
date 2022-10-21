## -----------------------------------------------------
# Function to simulate VAR model parameters
# for different dependence/correlation

sigphi <- function(p, rho)
{
  set.seed(1692, kind = "L'Ecuyer-CMRG" )
  A <- matrix(rnorm(p*p, mean = 0, sd = 1), p, p)
  B <- A %*% t(A)
  m <- max(eigen(B)$values)
  phi0 <- B/(m + 0.001)
  phik <- bdiag(rho*phi0)
  
  scratch <- diag((p)^2) - kronecker(phik, phik)
  V.s <- solve(scratch)%*%c(diag(p))
  V <- matrix(V.s, nrow = p, byrow = TRUE)
  Sigma <- solve(diag(p) - phik) %*% V + V %*% solve(diag(p) - phik) -V
  return(list(phi = phik, Sigma = Sigma))  
}


## -----------------------------------------------------
# Function to calculate the exact bias and variance

funBMexact <- function(b, x, y, r, c){
  n <- length(x)
  b <- floor(b)
  br <- floor(b/r)
  gamma0b1 <- 2*sum(x[2:b])
  gamma0b1r <- 2*sum(x[2:br])
  gamma0n1 <- 2*sum(x[2:n])
  gamma1b1 <- 2*sum((1:(b-1))*x[2:b])
  gamma1b1r <- 2*sum((1:(br-1))*x[2:br])
  gamma1n1 <- 2*sum(x[2:n]*seq(1, n-1))
  exact1 <- (1/(1-c))* (x[1] + (((n/b)*gamma0b1 - gamma0n1)/(n/b -1)) + 
            (gamma1n1 - (n/b)^2*gamma1b1)/(n*(n-b)/b)) - 
            (c/(1-c))* (x[1] + (((n/br)*gamma0b1r - gamma0n1)/(n/br -1)) +
            (gamma1n1 - (n/br)^2*gamma1b1r)/(n*(n-br)/br))
  bias.exact <- y - exact1
  variance.exact <- (8*(y^2)*b/n) +  (6*(y^2)*b/(r*n)) + (8*y^2*(r-1)/(r*(r - b/n)))*(b^2/n^2)
  mse.exact <- bias.exact^2 + variance.exact
  return(mse.exact)
}


## -----------------------------------------------------
# Function to calculate the proposed bias and variance for BM estimators

funBMi <- function(b, x, y, r, c){
  n <- length(x)
  b <- floor(b)
  br <- floor(b/r)
  gamma0b1 <- 2*sum(x[2:b])
  gamma0b1r <- 2*sum(x[2:br])
  gamma0n1 <- 2*sum(x[2:n])
  gamma1b1 <- 2*sum((1:(b-1))*x[2:b])
  gamma1b1r <- 2*sum((1:(br-1))*x[2:br])
  bias <- (1/(1-c))* ((n/(n-b))*(gamma0n1 - gamma0b1 + (gamma1b1)/b)) -
    (c/(1-c)* ((n/(n-br))*(gamma0n1 - gamma0b1r + (gamma1b1r)/br)))
  var <- (2*(y^2)*b/n) * (1/r + (r-1)/(r*(1-c)^2)) + 
    (2*(y^2)*(b/n)^2) * ((2*c)/((1-c)^2*r^2) - 2/r) 
  mse.bm <- bias^2 + var
  return(mse.bm)
}


## -----------------------------------------------------
# Function to calculate the proposed bias and variance for OBM estimators

funOBMi <- function(b, x, y, r, c){
  b <- floor(b)
  n <- length(x)
  br <- floor(b/r)
  gamma0b1 <- 2*sum(x[2:b])
  gamma0nb <- 2*sum(x[2:(n-b+1)])
  gamma0nbr <- 2*sum(x[2:(n-br+1)])
  gamma0b1r <- 2*sum(x[2:br])
  gamma0n1 <- 2*sum(x[2:n])
  gamma1b1 <- 2*sum((1:(b-1))*x[2:b])
  gamma1b1r <- 2*sum((1:(br-1))*x[2:br])
  bias <- (1/(1-c))*(gamma0n1 - gamma0b1 + (n*gamma1b1)/(b*(n-b)) - 
                       (b*n)/((n-b)*(n-b+1))*(gamma0b1 - gamma0nb)) - 
    (c/(1-c)*(gamma0n1 - gamma0b1r + (n*gamma1b1r)/(br*(n-br)) -
                (br*n)/((n-br)*(n-br+1))*(gamma0b1r - gamma0nbr)))
  var <- ((4/3)*(y^2)*b/n) * (1/r + (r-1)/(r*(1-c)^2)) + 
    ((4/3)*(y^2)*(b/n)^2) * ((2*c)/((1-c)^2*r^2) - 2/r) 
  mse.obm <- bias^2 + var
  return(mse.obm)
}


## -----------------------------------------------------
# current
funCurrbm <- function(b, x, y, r, c)
{
  mse.bm.curr <-  ((x/b)*(1-r*c)/(1-c))^2 + 
                  b/n*(y^2)*(1/r + (r-1)/(r*(1-c)^2)) 
                  # (8*(y^2)*b/n) +  (6*(y^2)*b/(r*n)) + 
                  # (8*y^2*(r-1)/(r*(r - b/n)))*(b^2/n^2)
  mse.bm.curr <- na.exclude(mse.bm.curr)
  return(mse.bm.curr)
}

funCurrobm <- function(b, x, y, r, c){
  mse.obm.curr <- ((x/b)*(1-r*c)/(1-c))^2 + 
                  ((4/3)*(y^2)*b/n) * (1/r + (r-1)/(r*(1-c)^2)) + 
                  ((4/3)*(y^2)*(b/n)^2) * ((2*c)/((1-c)^2*r^2) - 2/r) 
  mse.obm.curr <- na.exclude(mse.obm.curr)
  return(mse.obm.curr)
}

