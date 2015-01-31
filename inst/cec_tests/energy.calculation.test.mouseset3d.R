.testname <- "Energy calculation and covariances (mouseset3d)"

setup <- function()
{
  B <- as.matrix(read.table(system.file("cec_tests", "mouse3d.data", package="CEC")))
  C <- as.matrix(read.table(system.file("cec_tests", "centers3d.data", package="CEC")))
  C4 <- as.matrix(read.table(system.file("cec_tests", "centers43d.data", package="CEC")))
}

test.type.covariance <- function()
{
  given.cov = matrix(c(0.770118878, 0.005481129, -0.005991149, 0.005481129, 0.766972716, 0.008996509, -0.005991149, 0.008996509,  0.821481768), 3, 3)
  
  CE <- cec(B, centers=C, type="cov", param = given.cov, iter.max=20)
  
  cls = CE$nclusters[CE$iterations + 1]  
  expected.energy <- 0
  
  for(i in 1:cls)
  {
    cdata = CE$data[which(CE$cluster == i),] 
    cov <- cov.mle(cdata)
    Hx <- CEC:::H.covariance(cov, given.cov)
    p <- nrow(cdata) / nrow(B)    
    expected.energy <- expected.energy + p * (-log(p) + Hx)
  }   
  CEC:::checkNumericVectorEquals(expected.energy, CE$cost[CE$iterations + 1], msg="Energy")
}

test.type.fixedr.mixture <- function()
{
  r <- c(0.2, 0.3, 0.4) 
  
  CE <- cec(B, centers=C, type=c("fi", "fi", "fi"), param = r)
  
  cls = CE$nclusters[CE$iterations + 1]    
  expected.energy <- 0
  
  for(i in 1:cls)
  {
    cdata = CE$data[which(CE$cluster == i),] 
    cov <- cov.mle(cdata)
    Hx <- CEC:::H.fixedr(cov, r[i])
    p <- nrow(cdata) / nrow(B)    
    expected.energy <- expected.energy + p * (-log(p) + Hx)
  }  
  CEC:::checkNumericVectorEquals(expected.energy, CE$cost[CE$iterations + 1], msg="Energy")
}


test.type.spherical.one.cluster.removed <- function()
{  
  CE <- cec(B, C4, type="sp")
  
  cls = CE$nclusters[CE$iterations + 1]
  CEC:::checkNumericVectorEquals(3, cls, "Number of clusters")  
  expected.energy <- 0
  
  for(i in 1:cls)
  {
    cdata = CE$data[which(CE$cluster == i),] 
    cov <- cov.mle(cdata)
    Hx <- CEC:::H.spherical(cov)
    p <- nrow(cdata) / nrow(B)    
    expected.energy <- expected.energy + p * (-log(p) + Hx)
  }  
  CEC:::checkNumericVectorEquals(expected.energy, CE$cost[CE$iterations + 1], msg="Energy")
}

test.type.diagonal.spherical.mixture <- function()
{  
  CE <- cec(B, C, type=c("diag", "diag", "sp"))  
  
  cls = CE$nclusters[CE$iterations + 1]  
  CEC:::checkNumericVectorEquals(3, cls, "Number of clusters")  
  Hs <- list(CEC:::H.diagonal, CEC:::H.diagonal, CEC:::H.spherical) 
  
  expected.energy <- 0
  
  for(i in 1:cls)
  {
    cdata = CE$data[which(CE$cluster == i),] 
    cov <- cov.mle(cdata)
    Hx <- Hs[[i]](cov)
    p <- nrow(cdata) / nrow(B)    
    expected.energy <- expected.energy + p * (-log(p) + Hx)
  }  
  CEC:::checkNumericVectorEquals(expected.energy, CE$cost[CE$iterations + 1], msg="Energy")
}

test.type.eigenvalues.all.fixedr.mixture <- function()
{ 
  evals1 <- c(0.8240634, 0.7739987, 0.7595220)
  evals2 <- c(0.7240634, 0.5739987, 0.3595220)
  r <- 1.0
  
  evalsv = list(evals1, evals2)
  
  CE <- cec(B, C4, type=c("all", "eigen", "fixedr", "eigen"), param=list(evals1, r, evals2))
  
  cls = CE$final.nclusters
  CEC:::checkNumericVectorEquals(3, cls, "Number of clusters")
  Hs <- list(CEC:::H.all, CEC:::H.eigenvalues, CEC:::H.eigenvalues)  
  expected.energy <- 0
  
  for(i in 1:cls)
  {
    cdata = CE$data[which(CE$cluster == i),] 
    cov <- cov.mle(cdata)
    if (i == 1)
      Hx <- Hs[[i]](cov)
    else
      Hx <- Hs[[i]](cov, evalsv[[i-1]])
    
    p <- nrow(cdata) / nrow(B)    
    expected.energy <- expected.energy + p * (-log(p) + Hx)
  }
  CEC:::checkNumericVectorEquals(expected.energy, CE$cost[CE$iterations + 1], msg="Energy")
}
