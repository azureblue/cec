.testname <- "Energy calculation (ball1)"

setup <- function()
{
  B <- as.matrix(read.table(system.file("cec_tests", "ball1.data", package="CEC")))
  C <- as.matrix(read.table(system.file("cec_tests", "centers1.data", package="CEC")))
}


test.type.covariance <- function()
{
  cov <- CEC:::cov.mle(B)
  given.cov = matrix(c(2,1,1,3), 2,2)  
  
  expected.energy <- - log(1) + CEC:::H.covariance(cov, given.cov)
  
  CE <- cec(B, centers=1, type="cov", param = given.cov, iter.max=0)
  
  CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}

test.type.fixedr <- function()
{
  cov <- cov.mle(B)
  r <- 1.5
  
  expected.energy <- - log(1) + CEC:::H.fixedr(cov, r) 
  
  CE <- cec(B, centers=1, type="fix", param = 1.5, iter.max=0)
  
  CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}

test.type.spherical <- function()
{
  cov <- cov.mle(B)
  
  expected.energy <- - log(1) + CEC:::H.spherical(cov)
  
  CE <- cec(B, centers=1, type="sp", iter.max=0)
  
  CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}


test.type.diagonal <- function()
{
  cov <- cov.mle(B)

  expected.energy <- - log(1) + CEC:::H.diagonal(cov)
  
  CE <- cec(B, centers=1, type="diag", iter.max=0)
  
  CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}

test.type.eigenvalues <- function()
{
  cov <- cov.mle(B)
  evals <- c(0.1, 0.22)
  
  expected.energy <- CEC:::H.eigenvalues(cov, evals)
  
  CE <- cec(B, centers=1, type="eigen", param=evals, iter.max=0)
  
  CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}

test.type.eigenvalues <- function()
{
  cov <- cov.mle(B)
  evals <- c(0.1, 0.22)
  expected.energy <- - log(1) + CEC:::H.eigenvalues(cov, evals)
  
  CE <- cec(B, centers=1, type="eigen", param=evals, iter.max=0)
  
  CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}

test.type.all <- function()
{
  cov <- cov.mle(B)
  
  expected.energy <- CEC:::H.all(cov)
  
  CE <- cec(B, centers=1, type="all", iter.max=0)
  
  CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}

####################################################################################################################

test.type.spherical.cluster.removing <- function()
{
  cov <- cov.mle(B)
   
  expected.energy <- CEC:::H.spherical(cov)
  
  CE <- cec(B, C, type="sp", iter.max=20)
  
  CEC:::checkNumericVectorEquals(expected.energy, CE$cost[CE$iterations], msg="Energy")

}


