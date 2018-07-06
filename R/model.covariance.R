model.covariance <- function(type, cov, param)
{
    if (length(which(is.na(cov))) > 0)
        matrix(NA, nrow(cov), ncol(cov))
    else if (type == resolve.type("covariance"))
        param$cov
    else if (type == resolve.type("fixedr"))
        diag(ncol(cov)) * param$r
    else if (type == resolve.type("spherical"))
        diag(ncol(cov)) * sum(diag(ncol(cov)) * cov) / ncol(cov)
    else if (type == resolve.type("diagonal"))
        cov * diag(ncol(cov))
    else if (type == resolve.type("eigenvalues"))
    {    
        V <- eigen(cov)$vec
        D <- diag(sort(param$eigenvalues, decreasing=T))
        V %*% D %*% t(V)
    }
    else if (type == resolve.type("all"))
        cov
}