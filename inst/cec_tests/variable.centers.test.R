testname <- "Variable centers number"

setup <- function() {
    X1 = matrix(rnorm(1000), 500, 2)
    X2 = rbind(X1, X1 + 5)
    X3 = rbind(X2, X1 + 10)
}

test.should.use.1.cluster <- function() {
    C = cec(X1, 1:3, "sp", keep.removed=T, iter.max=0, card.min=4, nstart=10)
    CEC:::checkNumericEquals(1, nrow(C$centers))
}

test.should.use.2.cluster <- function() {
    C = cec(X2, 1:3, "sp", keep.removed=T, iter.max=0, card.min=4, nstart=10)
    CEC:::checkNumericEquals(2, nrow(C$centers))
}

test.should.use.3.cluster <- function() {
    C = cec(X3, 1:3, "sp", keep.removed=T, iter.max=0, card.min=4, nstart=5)
    CEC:::checkNumericEquals(3, nrow(C$centers))
}
