init.centers <- function(x, k, method = c("kmeans++", "random"))
{
    on.exit(.Call(release_cec_mem_r));
    method.int <- switch ( match.arg(method), "kmeans++" = 0, "random" = 1)
    
    if (! is.matrix(x)) stop("init.centers: x is not a matrix")
    if (k < 0) stop("init.centers: k < 0")
    
    n <- ncol(x)
    m <- nrow(x)

    centers <- .Call(cec_init_centers_r, x, k, method.int);
    centers  
}