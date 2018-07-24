testname <- "Split method"

setup <- function() {
    data("threeGaussians", package = "CEC")
    data("fourGaussians", package = "CEC")
    data("mixShapes", package = "CEC")

    mixShapesReduced = mixShapes[seq(1, nrow(mixShapes), 2),]
}

test.should.split.to.4.cluster <- function() {
    expected.cost = 2.530237
    tolerance = 0.001
    C = cec(fourGaussians, nstart = 5)
    CEC:::checkNumericEquals(4, C$nclusters)
    CEC:::checkNumericEquals(expected.cost, C$cost, msg = "cost", tolerance = tolerance)
}

test.should.split.to.7.cluster <- function() {
    expected.cost = 10.16551
    tolerance = 0.001
    C = cec(mixShapesReduced, 2, nstart = 2, split = T)
    CEC:::checkNumericEquals(7, C$nclusters)
    CEC:::checkNumericEquals(expected.cost, C$cost, msg = "cost", tolerance = tolerance)
}

test.should.limit.split.to.4.cluster <- function() {
    C = cec(mixShapesReduced, 1, nstart = 2, split = T, split.limit = 4)
    CEC:::checkNumericEquals(4, C$nclusters)
}

test.should.split.to.3.cluster.fixed.mean <- function() {
    C = cec(threeGaussians,, "mean", param = c(0, 0), nstart = 8)
    CEC:::checkNumericEquals(3, C$nclusters)
    CEC:::checkNumericEquals(1.726595, C$cost, tolerance = 0.00001)
}
