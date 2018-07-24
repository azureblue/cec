testname <- "Clustering"
setup <- function()
{ 
    data("fourGaussians", package = "CEC")
    centers <- as.matrix(read.table(system.file("cec_tests", "four.gaussians.centers.data", package="CEC")))
    expected <- dget(system.file("cec_tests", "four.gaussians.result.dp", package="CEC"))  
}

test.clustering.four.gaussians <- function()
{  
    CEC <- cec(fourGaussians, centers)
    CEC:::checkNumericVectorEquals(expected$cluster, CEC$cluster, msg="Clustering vector")
    CEC:::checkNumericVectorEquals(expected$cost, CEC$cost, msg="Energy")
    CEC:::checkNumericMatrixEquals(expected$centers, CEC$centers, msg="Centers")
    CEC:::checkNumericMatrixEquals(fourGaussians, CEC$data, msg="Data")
    CEC:::checkNumericVectorEquals(expected$probability, CEC$probability, msg="Probability")
    CEC:::checkNumericVectorEquals(expected$nclusters, CEC$nclusters, msg="Number of clusters")
    CEC:::checkNumericVectorEquals(expected$iterations, CEC$iterations, msg="Iterations")
}
