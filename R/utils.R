ball <- function(n = 4000, r = 1, dim = 2){ 
  M <- matrix(0,n,dim)
  count <- 0;
  rr <- r ^ 2
  while (count < n) 
  {
    p <- runif (dim, -r, r)        
    if (sum(p ^ 2) <= rr) M[count <- count + 1,] <- p      
  }
  M
}
