initcenters <- function(x, c, method = c("kmeans++", "random"))
{
  method.int <- switch ( match.arg(method), "kmeans++" = 0, "random" = 1)
  
  if (! is.matrix(x)) stop("initcenters: x is not a matrix")
  
  n <- ncol(x)
  m <- nrow(x)
  
  centers <- NULL
  
  if (method.int == 0)
    centers <- .Call(init_kmeanspp_r , x, c)
  else if (method.int == 1)
    centers <- .Call(init_random_r , x, c)
  
  centers  
}