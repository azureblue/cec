testname <- "Split method"

setup <- function() {
    data("fourGaussians", package = "CEC")
    data("mixShapes", package = "CEC")
}

test.should.split.to.4.cluster <- function() {
    expected.cost = 2.530237
    tolerance = 0.001
    C = cec(fourGaussians)
    CEC:::checkNumericEquals(4, C$nclusters)
    CEC:::checkNumericEquals(expected.cost, C$cost, msg = "cost", tolerance = tolerance)
}

test.should.split.to.7.cluster <- function() {
    expected.cost = 10.14958
    tolerance = 0.001
    C = cec(mixShapes, 2, split = T)
    CEC:::checkNumericEquals(7, C$nclusters)
    CEC:::checkNumericEquals(expected.cost, C$cost, msg = "cost", tolerance = tolerance)
}

test.should.limit.split.to.4.cluster <- function() {
    C = cec(mixShapes, 1, split = T, split.limit = 4)
    CEC:::checkNumericEquals(4, C$nclusters)
}
