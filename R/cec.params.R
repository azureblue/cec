resolve.type <- function(type)
{
  types <- c("covariance", "fixedr", "spherical", "diagonal", "all")
  int.type <- switch (match.arg(type, types), covariance = 0, fixedr = 1, spherical = 2, diagonal = 3, all = 4)
  int.type  
}

create.cec.params <- function(k, n, type, param)
{
  params <- rep(list(NA), k)
  
  if (length(type) == 1) 
  {      
    type.i <- resolve.type(type)
    type.v  <-  rep(resolve.type(type), k)
    if (type.i == resolve.type("covariance")) 
    {
      if (!hasArg(param)) stop("Missing requierd argument: 'param'.")
      
      if (is.matrix(param))
        cov <- param
      else if (is.list(param))
        cov <- param[[1]]
      
      if (!is.array(cov)) stop("Illegal argument: illegal parameter for \"covariance\" type.")    
      if (ncol(cov) != n) stop("Illegal argument: illegal parameter for \"covariance\" type.")    
      if (nrow(cov) != n) stop("Illegal argument: illegal parameter for \"covariance\" type.")
      if("try-error" %in% class(try(chol(cov), silent=TRUE))) 
        stop("Illegal argument: illegal parameter for \"covariance\" type. Matrix must be positive-definite.")
      i_cov = solve(cov)  
      params <- rep(list(list(cov, i_cov)), k)
    }
    else if ( type.i == resolve.type("fixedr")) 
    {
      r <- param[1]
      if (!r > 0)  stop("Illegal argument: illegal parameter for \"fixedr\" type.")
      params = rep(list(r), k)
    }   
  }
  else 
  {
  if (k != length(type)) 
    stop("Illegal argument: illegal length of \"type\" vector.")
  
  if (hasArg(param)) 
    param <- param[!param %in% list(NULL, NA)]
  
  idx <- 0
     for (i in 1:length(type))
     {
       type.i <- resolve.type(type[i])
       if (type.i == resolve.type("covariance")) 
       {
         idx <- idx + 1
         
         #check param length
         if (length(param) < idx)
           stop("Illegal argument: illegal param length.")

         #expect matrix param
         cov <- param[[idx]]
         
         if (!is.array(cov)) stop("Illegal argument: illegal parameter for \"covariance\" type.")    
         if (ncol(cov) != n) stop("Illegal argument: illegal parameter for \"covariance\" type.")    
         if (nrow(cov) != n) stop("Illegal argument: illegal parameter for \"covariance\" type.")
         if("try-error" %in% class(try(chol(cov), silent=TRUE))) 
           stop("Illegal argument: illegal parameter for \"covariance\" type. Matrix must be positive-definite.")
         i.cov = solve(cov)  
       params[[i]] <- list(cov, i.cov)
       }
       else if (type.i == resolve.type("fixed")) 
       {
         idx <- idx + 1
         
         #check param length
         if (length(param) < idx)
           stop("Illegal argument: illegal param length.")
         
         #expect numeric param
         r = param[[idx]]
         if (!is.numeric(r)) stop("Illegal argument: illegal parameter for \"fixedr\" type.")
         if (!r > 0)  stop("Illegal argument: illegal parameter for \"fixedr\" type.")
         params[i] <- r
       }
     }          
  }
  params
}