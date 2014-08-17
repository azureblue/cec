.testname <- "Energy calculation (ball1)"

setup <- function()
{
  B <- as.matrix(read.table(system.file("cec_tests", "ball1.data", package="cec")))
}


test.type.covariance <- function()
{
  bcov <- cov.mle(B)
  given.cov = matrix(c(2,1,1,3), 2,2)
  igiven.cov <- solve(given.cov)
  Hx <- 2 / 2 * log(2 * pi) + 1/2 * sum(diag(igiven.cov %*% bcov)) + 1 / 2 * log (det(given.cov))  
  expected.energy <- - log(1) + Hx
  
  CE <- cec(B, centers=1, type="cov", param = given.cov, iter.max=0)
  
  checkNumericVectorEquals(expected.energy, CE$energy[1], msg="Energy")
}

test.type.fixedr <- function()
{
  bcov <- cov.mle(B)
  r <- 1.5
  Hx <- log(2 * pi) + 1/(2*r) * sum(diag(bcov)) + log (r)  
  expected.energy <- Hx  
  
  CE <- cec(B, centers=1, type="fix", param = 1.5, iter.max=0)
  
  checkNumericVectorEquals(expected.energy, CE$energy[1], msg="Energy")
}

test.type.spherical <- function()
{
  bcov <- cov.mle(B)
  Hx <- log(2 * pi * 2.718281828 / 2) + log (sum(diag(bcov)))
  expected.energy <- Hx  
  
  CE <- cec(B, centers=1, type="sp", iter.max=0)
  
  checkNumericVectorEquals(expected.energy, CE$energy[1], msg="Energy")
}


test.type.diagonal <- function()
{
  bcov <- cov.mle(B)
  Hx <- log(2 * pi * 2.718281828) + log(bcov[1,1] * bcov[2,2]) / 2
  expected.energy <- Hx  
  
  CE <- cec(B, centers=1, type="diag", iter.max=0)
  
  checkNumericVectorEquals(expected.energy, CE$energy[1], msg="Energy")
}

test.type.all <- function()
{
  bcov <- cov.mle(B)
  Hx <- log(2 * pi * 2.718281828) + log(det(bcov)) / 2
  expected.energy <- Hx  
  
  CE <- cec(B, centers=1, type="all", iter.max=0)
  
  checkNumericVectorEquals(expected.energy, CE$energy[1], msg="Energy")
}