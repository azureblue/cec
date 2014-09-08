
cec <- function(
  x, 
  centers,   
  iter.max      = 20,
  nstart        = 1,
  centers.init  = c("kmeans++", "random"), 
  type          = c("covariance", "fixedr", "spherical", "diagonal", "eigenvalues", "all"),  
  param,
  card.min      = "5%",
  keep.removed  = F,
  interactive   = F,
  readline      = T
  )
{ 
  ## checking arguments

  if (!hasArg(x)) stop("Missing requierd argument: 'x'.")
  if (!hasArg(centers)) stop("Missing requierd argument: 'centers'.")
  
  if (!is.matrix(x)) stop("Illegal argument: 'x' is not a matrix.")
  if (ncol(x) < 1) stop("Illegal argument: ncol(x) < 1.")
  if (nrow(x) < 1) stop("Illegal argument: nrow(x) < 1.")  
  
  if (!is.matrix(centers))
  {
    if (centers < 1) stop("Illegal argument: 'centers' < 1")
    centers.initilized <- F
  } 
  else 
  {
    if (ncol(x) != ncol(centers)) stop("Illegal argument: number of columns of 'x' and 'centers' are not equal.")
    if (nrow(centers) < 1) stop("Illegal argument: nrow(centers) < 1.")
    centers.initilized <- T
  }
  
  if (!(attr(regexpr("[\\.0-9]+%{0,1}", perl=TRUE, text=card.min), "match.length") == nchar(card.min)))
    stop("Illegal argument: 'card.min' in wrong format.")  
  
  
  if (hasArg(centers.init)) match.arg(centers.init)
  if (!hasArg(type))
    type <- "all"
  
  ####################################################
  
  # run interactive mode if requested
  
  if (interactive)
    return(cec_interactive(x, centers, iter.max, 1, centers.init, type, param, card.min, keep.removed, readline))
    
  ####################################################
  
  n = ncol(x)
  m = nrow(x)
  k <- NULL
  if (is.matrix(centers)) k <- nrow(centers) 
  else k <- centers

  params <- create.cec.params(k, n, type, param)
  type <- as.integer(vapply(type, resolve.type, 0))
  if (length(type) == 1) type <- rep(type, k)
    
  
  if (substr(card.min, nchar(card.min), nchar(card.min)) == "%") 
  {
    card.min = as.integer(as.double(substr(card.min , 1, nchar(card.min) - 1)) * m / 100)
  } 
  else
  {
    card.min = as.integer(card.min)
  }
  
    tenergy = .Machine$integer.max
    Z <- NULL
    if (!is.matrix(centers)) centers <- initcenters(x, centers, centers.init)        
    ok.flag <- F
    startTime <- proc.time()     
    for (start in 1:nstart) {
      if (!centers.initilized) centers <- initcenters(x, k, centers.init)     
      tryCatch( 
       {
         X <- .Call(cec_r, x, centers, iter.max, type, card.min, params)  
         ok.flag <- T
         if (X$iterations < 0 || X$energy[X$iterations + 1] < tenergy)
         {
           tenergy = X$energy[X$iterations+1]
           Z <- X
         }
       }, error = function(er) {
         warning(paste("Error at start #", start, ": ", er$message), immediate.=T, call.=F)
       }
     )     
    }
     if (ok.flag == F) 
     {
       stop("All starts faild with error.")
     }    
    execution.time = as.vector((proc.time() - startTime))[1]
  
    Z$centers[is.nan(Z$centers)] <- NA
  
    tab <- tabulate(Z$cluster)
    probability <- vapply(tab, function(c.card){c.card / m}, 0)
  
    if (!keep.removed)
    {
      cluster.map = 1:k
      na.rows = which(is.na(Z$centers[, 1]))
      if (length(na.rows) > 0)
      {
        for (i in 1:length(na.rows))        
          for (j in na.rows[i]:k)          
            cluster.map[j] <- cluster.map[j] - 1    
        
        Z$cluster  <- as.integer(vapply(Z$cluster,function(asgn) {as.integer(cluster.map[asgn])}, 0))
        
        #in case of having single row - convert it to a matrix
        Z$centers <- matrix(Z$centers[-na.rows,],,n)
        
        Z$covariances <- Z$covariances[-na.rows]
        probability <- probability[-na.rows]
        params <- params[-na.rows]
        type <- type[-na.rows]
      }     
    }
    structure(list(
      data              = x,
      cluster           = Z$cluster,
      centers           = Z$centers, 
      probability       = probability,
      energy            = Z$energy  [1:(Z$iterations+1)], 
      nclusters         = Z$nclusters,
      iterations        = Z$iterations, 
      time              = execution.time,
      covariances       = Z$covariances,
      covariances.model = lapply(list.triple(type, Z$covariances, params), model.covariance)
    ), class = "cec");    
}

list.triple <- function(type.v, cov.l, param.l)
{
  len <- length(type.v)
  L <- rep(list(NA), len)
  for (i in 1:len)
  {
    L[[i]] <- list(type.v[i], cov.l[[i]], param.l[[i]])
  }
  L
}

model.covariance <- function(type.cov.param)
{
  type  <- type.cov.param[[1]]
  cov   <- type.cov.param[[2]]
  param <- type.cov.param[[3]]
  
  if (length(which(is.na(cov))) > 0)
  {
    matrix(NA, nrow(cov), ncol(cov))
  }  
  else if (type == resolve.type("covariance"))
  {
    param[[1]]
  }  
  else if (type == resolve.type("fixedr"))
  {
    diag(ncol(cov)) * param
  }
  else if (type == resolve.type("spherical"))
  {
    diag(ncol(cov)) * sum(diag(ncol(cov)) * cov) / ncol(cov)
  }
  else if (type == resolve.type("diagonal"))
  {
    cov * diag(ncol(cov))
  }
  else if (type == resolve.type("eigenvalues"))
  {
    
    V <- eigen(cov)$vec
    D <- diag(sort(param, decreasing=T))
    
    V %*% D %*% t(V)
  }
  else if (type == resolve.type("all"))
  {
    cov
  }
}

print.cec <- function(x, ...)
{
  cat("CEC clustering result: \n\n")
  cat("Clustering vector: \n")
  print(x$cluster)
  cat("\nCluster means:\n")
  print(x$centers)
  cat("\nProbability vector:\n")
  print(x$probability)
  cat("\nEnergy at each iteration:\n")
  print(x$energy)  
  cat("\nNumber of clusters at each iteration:\n")
  print(x$nclusters)  
  cat("\nNumber of iterations:\n")
  print(x$iterations)
  cat("\nComputation time:\n")
  print(x$time)
  cat("\nAvailable components:\n")
  print(c("data", "cluster", "centers", "probabilities", "energy", "nclusters", "iterations", "covariances", "covariances.model", "time" ))
}

plot.cec <- function(x, col, cex = 0.5, pch = 16, cex.centers = 1, pch.centers = 8, ellipses.lwd = 4, ellipses = TRUE, model = FALSE, xlab = "x", ylab= "y", ...)
{
  if (ncol (x $ data) != 2 )
    stop("plotting available only for 2-dimensional data")
  
  if(!hasArg(col)) col = x$cluster;
  plot(x$data, col=col, cex = cex, pch = pch,  xlab = xlab, ylab = ylab, ...)    
  points(x$centers, cex = cex.centers, pch = pch.centers) 
  if (x$iterations > -1)
  if (ellipses)
  {    
    for (i in 1:nrow(x$centers))     
      if (! is.na(x$centers[i, 1]))
      {         
        err = FALSE        
        tryCatch(
          {
            cov <- NA
            if (model == T)
              cov <- x$covariances.model[[i]]
            else
              cov <- x$covariances[[i]]
            pts <- ellipse(x$centers[i, ], cov)
            lines(pts, lwd = ellipses.lwd)
          },
          warning = function(e) {warning("some ellipses will not be drawn (probably not positive-definite covariance matrix)")},
          error = function(e) {warning("some ellipses will not be drawn (probably not positive-definite covariance matrix)")}, 
          finally = {})     
      }
    }  
}

cec.plot.energy <- function(C, xlab="Iteration", ylab="Energy", lwd=5, col="red", lwd.points=5, pch.points=19, col.points="black", mgp=c(1.5,0.5,0), ...)
{
  plot(x = 1:(length(C$energy) - 1), y = C$energy[2:(length(C$energy))], xlab=xlab, ylab=ylab, type="l", lwd=lwd, col=col, mgp=mgp, ...)
  points(x = 1:(length(C$energy) - 1), y = C$energy[2:(length(C$energy))], lwd=lwd.points, pch=pch.points, col=col.points)
  title("Energy at each iteraion")
}


cec_interactive <- function(
  x, 
  centers,   
  iter.max      = 20,
  nstart        = 1,
  centers.init  = c("kmeans++", "random"), 
  type          = c("covariance", "fixedr", "spherical", "diagonal", "eigenvalues", "all"),  
  param,
  card.min      = "5%",
  keep.removed  = F,
  readline      = T
)
{
  par 
  {
    old.ask = par()["ask"]   
    n = ncol(x)
    if (n != 2) 
      stop("interactive mode available only for 2-dimensional data")    
    i <- -1    
    if (!is.matrix(centers)) centers <- initcenters(x, centers, centers.init)
    if (readline)
    {
      ignore = readline(prompt="After each iteration you may:\n - press <Enter> for next iteration \n - write number <n> (may be negative one) and press <Enter> for next <n> iterations \n - write 'q' and abort execution.\n Press <Return>.\n")        
      par(ask = FALSE)
    }
    else
    {
      par(ask = TRUE)
    }
    while (TRUE) 
    {      
      Z <- cec(x, centers, i, 1, centers.init, type, param, card.min, keep.removed, F);
      
      if(i > Z$iterations | i>= iter.max) 
        break
      
      desc = ""
      if (i == -1)
        desc = "(initial cluster centers and first assignment to groups)"
      else if (i == 0)
        desc = "(position of center means before first iteration)"      
      
      cat("Iterations:", Z$iterations, desc, "energy:", Z$energy[(Z$iterations + 1)]," \n ")
      
      if (i == -1)
        plot(Z, ellipses = FALSE)
      else
        plot(Z, ellipses = TRUE)
      
      if (readline) 
      {
        line = readline(prompt="Press <Enter> OR write number OR write 'q':");     
        lineint = suppressWarnings(as.integer(line))
        
        if (!is.na(lineint)) 
        {
          i = i + lineint - 1
          if (i < -1) 
            i = -2        
        } 
        else if (line == "q" | line == "quit") {
          break
        }
      }
      i = i + 1
    }    
    plot(Z, ellipses="TRUE")
    par(ask = old.ask)
    if (readline) 
    {
      ignore = readline(prompt="Press <Enter>:")
    }
    Z
  }
}
