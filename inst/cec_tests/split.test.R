testname <- "Split method"

setup <- function() {
    data("fourGaussians", package = "CEC")
    data("mixShapes", package = "CEC")

    mixShapesReduced = mixShapes[seq(1, nrow(mixShapes), 2),]
}

test.should.split.to.4.cluster <- function() {
    expected.cost = 2.530237
    tolerance = 0.001
    C = cec(fourGaussians, nstart=2)
    CEC:::checkNumericEquals(4, C$nclusters)
    CEC:::checkNumericEquals(expected.cost, C$cost, msg = "cost", tolerance = tolerance)
}

test.should.split.to.7.cluster <- function() {
    expected.cost = 10.16551
    tolerance = 0.001
    C = cec(mixShapesReduced, 2, split = T)
    CEC:::checkNumericEquals(7, C$nclusters)
    CEC:::checkNumericEquals(expected.cost, C$cost, msg = "cost", tolerance = tolerance)
}

test.should.limit.split.to.4.cluster <- function() {
    C = cec(mixShapesReduced, 1, split = T, split.limit = 4)
    CEC:::checkNumericEquals(4, C$nclusters)
}
