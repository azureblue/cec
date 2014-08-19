#source("tests")
.testname <- "Clustering"
setup <- function()
{
  M <- as.matrix(read.table(system.file("cec_tests", "mouse1.data", package="cec")))
  C <- as.matrix(read.table(system.file("cec_tests", "centers1.data", package="cec")))
  expected <- dget(system.file("cec_tests", "cec1.dp", package="cec"))
  CEC <- cec(M, C)
}
test.clustering.mouse1 <- function()
{  
  checkNumericVectorEquals(expected$cluster, CEC$cluster, msg="Clustering vector")
  checkNumericVectorEquals(expected$energy, CEC$energy, msg="Energy")
  checkNumericMatrixEquals(expected$centers, CEC$centers, msg="Centers")  
}