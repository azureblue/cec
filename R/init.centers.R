init.centers <- function(x, k, method = c("kmeans++", "random"))
{
    on.exit(.Call(release_cec_mem_r));
    method <- switch ( match.arg(method), "kmeans++" = "kmeanspp", "random" = "random")
    
    if (! is.matrix(x)) stop("init.centers: x is not a matrix")
    if (k < 0) stop("init.centers: k < 0")

    centers <- .Call(cec_init_centers_r, x, k, method);
    centers  
}