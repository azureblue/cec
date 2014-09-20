.testname <- "Energy calculation (mouseset1)"

setup <- function()
{
  B <- as.matrix(read.table(system.file("cec_tests", "mouse1.data", package="CEC")))
}

test.type.covariance <- function()
{
  cov <- cov.mle(B)
  given.cov = matrix(c(2,1,1,3), 2,2)
  igiven.cov <- solve(given.cov)
  
  expected.energy <- CEC:::H.covariance(cov, given.cov)  
  
  CE <- cec(B, centers=1, type="cov", param = given.cov, iter.max=0)
  
  CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}

test.type.fixedr <- function()
{
  cov <- cov.mle(B)
  r <- 1.5
  
  expected.energy <- CEC:::H.fixedr(cov, r)
  
  CE <- cec(B, centers=1, type="fix", param = 1.5, iter.max=0)
  
  CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}

test.type.spherical <- function()
{
  bcov <- cov.mle(B)
  Hx <- log(2 * pi * 2.718281828 / 2) + log (sum(diag(bcov)))
  expected.energy <- Hx  
  
  CE <- cec(B, centers=1, type="sp", iter.max=0)
  
  CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}


test.type.diagonal <- function()
{
  cov <- cov.mle(B)  
  
  expected.energy <- CEC:::H.diagonal(cov)
  
  CE <- cec(B, centers=1, type="diag", iter.max=0)
  
  CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}

test.type.all <- function()
{
  cov <- cov.mle(B)
  
  expected.energy <- CEC:::H.all(cov)
  
  CE <- cec(B, centers=1, type="all", iter.max=0)
  
  CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}
