
cec <- function(
  x, 
  centers,   
  type          = c("covariance", "fixedr", "spherical", "diagonal", "eigenvalues", "all"),  
  iter.max      = 25,
  nstart        = 1,
  param,
  centers.init  = c("kmeans++", "random"), 
  card.min      = "5%",
  keep.removed  = F,
  interactive   = F,
  readline      = T
)
{ 
  # check arguments  
  if (!hasArg(x)) stop("Missing requierd argument: 'x'.")
  if (!hasArg(centers)) stop("Missing requierd argument: 'centers'.")
  
  if (!is.matrix(x)) stop("Illegal argument: 'x' is not a matrix.")
  if (ncol(x) < 1) stop("Illegal argument: ncol(x) < 1.")
  if (nrow(x) < 1) stop("Illegal argument: nrow(x) < 1.")

  if (!all(complete.cases(x))) stop("Illegal argument: 'x' contains NA values.")
  if (!all(complete.cases(centers))) stop("Illegal argument: 'centers' contains NA values.")
  
  var.centers <- F
  if (!is.matrix(centers))
  {
    if (length(centers) > 1)
    {
      var.centers <- T
      if (hasArg(nstart) && nstart != length(centers))
        stop("Illegal argument: 'nstart' != length(centers)")
      nstart <- length(centers)
    }
    for (i in centers) 
      if (i < 1) stop("Illegal argument: 'centers' < 1")    
    centers.initilized <- F
  } 
  else 
  {
    if (ncol(x) != ncol(centers)) stop("Illegal argument: ncol(x) != ncol(centers).")
    if (nrow(centers) < 1) stop("Illegal argument: nrow(centers) < 1.")
    centers.initilized <- T
  }
  
  if (!(attr(regexpr("[\\.0-9]+%{0,1}", perl=TRUE, text=card.min), "match.length") == nchar(card.min)))
    stop("Illegal argument: 'card.min' in wrong format.")  
  
  
  if (hasArg(centers.init)) match.arg(centers.init)
  if (!hasArg(type))
    type <- "all"
  
  if (var.centers && length(type) > 1)
    stop("Illegal argument: 'type' with length > 1 is not supported for variable number of centers ('centers' as a vector).")
  
  # run interactive mode if requested  
  if (interactive)
    return(cec_interactive(x, centers, type, iter.max, 1, param, centers.init, card.min, keep.removed, readline))
  
  n = ncol(x)
  m = nrow(x)
  k <- NULL
  if (is.matrix(centers)) 
    k <- nrow(centers) 
  else 
    k <- centers
  
  if (substr(card.min, nchar(card.min), nchar(card.min)) == "%") 
    card.min = as.integer(as.double(substr(card.min , 1, nchar(card.min) - 1)) * m / 100)
  else
    card.min = as.integer(card.min)
  
  # card.min must be greater than the dimension of the data
  card.min = max(card.min, n + 1)
  
  tenergy = .Machine$integer.max
  Z <- NULL
  ok.flag <- F    
  
  # prepare input for C function cec_r
  if (!var.centers)
  {
    params <- create.cec.params(k, n, type, param)
    types <- as.integer(vapply(type, resolve.type, 0))
    if (length(types) == 1) 
      types <- rep(types, k)
  }
  
  startTime <- proc.time()     
  
  # main loop
  for (start in 1:nstart) 
  {      
    if (var.centers)
    {
      # prepare input for C function cec_r with regards to the variable number of centers
      k <- centers[start]        
      params <- create.cec.params(k, n, type, param)
      types <- as.integer(vapply(type, resolve.type, 0))
      if (length(types) == 1) 
        types <- rep(types, k)
    }
    # generate initial centers or use provided centers matrix
    if (!centers.initilized) 
      centers.matrix <- initcenters(x, k, centers.init)     
    else
      centers.matrix <- centers      
    
    tryCatch(
      {
        # perform the clustering by calling C function cec_r
        X <- .Call(cec_r, x, centers.matrix, iter.max, types, card.min, params)  
        ok.flag <- T
        if (X$iterations < 0 || X$energy[X$iterations + 1] < tenergy)
        {
          # keep the clustering results with the lowest energy (cost function)
          tenergy = X$energy[X$iterations+1]
          Z <- X
          k.final <- k
          types.final <- types
          params.final <- params
        }
      }, error = function(er) {
        warning(paste("Error at start #", start, ": ", er$message), immediate.=T, call.=F)
      })     
  }
  
  if (ok.flag == F) 
  {
    stop("All starts faild with error.")
  }    
   
  # prepare the results  
  
  execution.time = as.vector((proc.time() - startTime))[1]
  
  Z$centers[is.nan(Z$centers)] <- NA
  
  tab <- tabulate(Z$cluster)
  probability <- vapply(tab, function(c.card){c.card / m}, 0)
  
  # change cluster assignment if keep.removed == F
  if (!keep.removed)
  {
    cluster.map = 1:k.final
    na.rows = which(is.na(Z$centers[, 1]))
    if (length(na.rows) > 0)
    {
      for (i in 1:length(na.rows))        
        for (j in na.rows[i]:k.final)          
          cluster.map[j] <- cluster.map[j] - 1    
      
      Z$cluster  <- as.integer(vapply(Z$cluster,function(asgn) {as.integer(cluster.map[asgn])}, 0))
      
      # in case of having single row - convert it to a matrix
      Z$centers <- matrix(Z$centers[-na.rows,],,n)        
      Z$covariances <- Z$covariances[-na.rows]
      probability <- probability[-na.rows]
      params.final <- params.final[-na.rows]
      types.final <- types.final[-na.rows]
    }     
  }
  covs = length(types.final)
  covariances.model = rep(list(NA), covs)
  
  # obtain the covariances of the model
  for(i in 1:covs)
    covariances.model[[i]] = model.covariance(types.final[[i]], Z$covariances[[i]], params.final[[i]])

  structure(list(
    data                = x,
    cluster             = Z$cluster,
    centers             = Z$centers, 
    probability         = probability,
    cost.function       = Z$energy[1:(Z$iterations + 1)], 
    nclusters           = Z$nclusters,
    final.cost.function = Z$energy[Z$iterations + 1], 
    final.nclusters     = Z$nclusters[Z$iterations + 1],
    iterations          = Z$iterations, 
    time                = execution.time,
    covariances         = Z$covariances,
    covariances.model   = covariances.model
  ), class = "cec");    
}

model.covariance <- function(type, cov, param)
{
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
  cat("\nProbability vector:\n")
  print(x$probability)
  cat("\nMeans of clusters:\n")
  print(x$centers)
  cat("\nCost function at each iteration:\n")
  print(x$cost)  
  cat("\nNumber of clusters at each iteration:\n")
  print(x$nclusters)  
  cat("\nNumber of iterations:\n")
  print(x$iterations)
  cat("\nComputation time:\n")
  print(x$time)
  cat("\nAvailable components:\n")
  print(c("data", "cluster", "probabilities", "centers", "cost.function", "nclusters", "final.cost.function", "final.nclusters", "iterations", "covariances", "covariances.model", "time" ))
}

plot.cec <- function(x, col, cex = 0.5, pch = 16, cex.centers = 1, pch.centers = 8, ellipses.lwd = 4, ellipses = TRUE, model = T, xlab = "x", ylab= "y", ...)
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
            #   warning = function(e) {warning("some ellipses will not be drawn (probably not positive-definite covariance matrix)")},
            #   error = function(e) {warning("some ellipses will not be drawn (probably not positive-definite covariance matrix)")}, 
            finally = {})     
        }
    }  
}

cec.plot.cost.function <- function(C, xlab="Iteration", ylab="Cost function", lwd=5, col="red", lwd.points=5, pch.points=19, col.points="black", mgp=c(1.5,0.5,0), ...)
{
  plot(x = 1:(length(C$cost) - 1), y = C$cost[2:(length(C$cost))], xlab=xlab, ylab=ylab, type="l", lwd=lwd, col=col, mgp=mgp, ...)
  points(x = 1:(length(C$cost) - 1), y = C$cost[2:(length(C$cost))], lwd=lwd.points, pch=pch.points, col=col.points)
  title("Cost function at each iteraion")
}

cec_interactive <- function(
  x, 
  centers,   
  type          = c("covariance", "fixedr", "spherical", "diagonal", "eigenvalues", "all"),  
  iter.max      = 40,
  nstart        = 1,
  param,
  centers.init  = c("kmeans++", "random"), 
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
    Z <- cec(x, centers, type, i, 1, param, centers.init , card.min, keep.removed, F);
    
    if(i > Z$iterations | i>= iter.max) 
      break
    
    desc = ""
    if (i == -1)
      desc = "(initial cluster centers and first assignment to groups)"
    else if (i == 0)
      desc = "(position of center means before first iteration)"      
    
    cat("Iterations:", Z$iterations, desc, "cost function:", Z$cost[(Z$iterations + 1)]," \n ")
    
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
