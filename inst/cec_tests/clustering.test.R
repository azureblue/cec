#source("tests")
.testname <- "Clustering"
setup <- function()
{
  M <- as.matrix(read.table(system.file("cec_tests", "mouse1.data", package="CEC")))
  C <- as.matrix(read.table(system.file("cec_tests", "centers1.data", package="CEC")))
  expected <- dget(system.file("cec_tests", "cec1.dp", package="CEC"))  
}
test.clustering.mouse1 <- function()
{  
  CE <- cec(M, C)
  CEC:::checkNumericVectorEquals(expected$cluster, CE$cluster, msg="Clustering vector")
  CEC:::checkNumericVectorEquals(expected$cost, CE$cost, msg="Energy")
  CEC:::checkNumericMatrixEquals(expected$centers, CE$centers, msg="Centers")  
}
