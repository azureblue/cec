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
    threads       = 1,
    split         = F,
    split.depth   = 8,
    split.tries   = 5,
    split.limit   = 100,
    split.initial.starts = 1,
    readline      = T
)
{
    # check arguments  
    if (!hasArg(x)) stop("Missing requierd argument: 'x'.")

    if (!hasArg(centers)) {
        centers <- 1
        split = T
    }
    
    if (iter.max < 0) stop("Illegal argument: iter.max < 0.")
    if (!is.matrix(x)) stop("Illegal argument: 'x' is not a matrix.")
    if (ncol(x) < 1) stop("Illegal argument: ncol(x) < 1.")
    if (nrow(x) < 1) stop("Illegal argument: nrow(x) < 1.")
    
    if (!all(complete.cases(x))) stop("Illegal argument: 'x' contains NA values.")
    if (!all(complete.cases(centers))) stop("Illegal argument: 'centers' contains NA values.")
    
    var.centers = NULL
    centers.mat = NULL
    
    if (!is.matrix(centers))
    {
        if (length(centers) > 1)
            var.centers <- centers
        else
            var.centers <- c(centers)
        
        for (i in centers) 
            if (i < 1) stop("Illegal argument: 'centers' < 1")

        centers.initilized <- F
    } 
    else 
    {
        if (ncol(x) != ncol(centers)) stop("Illegal argument: ncol(x) != ncol(centers).")
        if (nrow(centers) < 1) stop("Illegal argument: nrow(centers) < 1.")
        var.centers <- c(nrow(centers))
        centers.mat = centers
        centers.initilized <- T
    }
    
    if (!(attr(regexpr("[\\.0-9]+%{0,1}", perl=TRUE, text=card.min), "match.length") == nchar(card.min)))
        stop("Illegal argument: 'card.min' in wrong format.")  
    
    if (centers.initilized)
        init.method.name = "none"
    else if (hasArg(centers.init))
        init.method.name = switch (match.arg(centers.init), "kmeans++" = "kmeanspp", "random" = "random")
    else init.method.name = "kmeanspp"
    
    
    if (!hasArg(type))
        type <- "all"
    
    if (length(type) > 1 && length(var.centers) > length(type))
        stop("Illegal argument: 'type' with length > 1 should be equal or greater than the length of the vector of variable number of centers ('centers' as a vector).")
    
    # run interactive mode if requested  
    if (interactive) {
        if (split == T) stop("Interactive mode not available in split mode")
        if (length(var.centers) > 1) stop("Interactive mode not available for variable centers")
        if (nstart > 1) stop("Interactive mode not available for multiple starts")
        return(cec.interactive(x, centers, type, iter.max, 1, param, centers.init, card.min, keep.removed, readline))
    }
    
    n = ncol(x)
    m = nrow(x)
    
    if (substr(card.min, nchar(card.min), nchar(card.min)) == "%") 
        card.min = as.integer(as.double(substr(card.min , 1, nchar(card.min) - 1)) * m / 100)
    else
        card.min = as.integer(card.min)
    
    # card.min must be greater than the dimension of the data
    card.min = max(card.min, n + 1)

    # prepare input for C function cec_r
    k <- max(var.centers)
    
    startTime <- proc.time()     
    
    centers.r = list(
        init.method = init.method.name,
        var.centers = as.integer(var.centers),
        mat = centers.mat
    )
    
    if (threads == "auto")
        threads = 0
    
    control.r = list(
        min.card = as.integer(card.min),
        max.iters = as.integer(iter.max),
        starts = as.integer(nstart),
        threads = as.integer(threads)
    )
    
    models.r = create.cec.params.for.models(k, n, type, param)
    
    # perform the clustering by calling C function cec_r
    if (split) {
        for (i in 1:k)
            if (models.r[[i]]$type != models.r[[1]]$type)
                stop("mixing model types is currently not supported in split mode")
        
        models.r
        split.r = list(
            depth = as.integer(split.depth),
            limit = as.integer(split.limit),
            tries = as.integer(split.tries),
            initial.starts = as.integer(split.initial.starts)
        )
        Z <- .Call(cec_split_r, x, centers.r, control.r, models.r, split.r)
    }
    else
        Z <- .Call(cec_r, x, centers.r, control.r, models.r)
    
    k.final <- nrow(Z$centers)
    
    # prepare the results  
    
    execution.time = as.vector((proc.time() - startTime))[3]
    
    Z$centers[is.nan(Z$centers)] <- NA
    
    tab <- tabulate(Z$cluster)
    probability <- vapply(tab, function(c.card){c.card / m}, 0)
    
    #TODO: change this temporary hack
    model.one = models.r[[1]];
    
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
                models.r <- models.r[-na.rows]
        }     
    }
    covs = length(Z$covariances)
    covariances.model = rep(list(NA), covs)
    
    #TODO: change this temporary hack
    if (split)
        models.r = rep(list(model.one), covs)
    
    # obtain the covariances of the model
    for(i in 1:covs)
        covariances.model[[i]] = model.covariance(models.r[[i]]$type, Z$covariances[[i]], models.r[[i]]$params)
    
    structure(list(
        data                = x,
        cluster             = Z$cluster,
        centers             = Z$centers, 
        probability         = probability,
        cost.function       = Z$energy, 
        nclusters           = Z$nclusters,
        iterations          = Z$iterations, 
        time                = execution.time,
        covariances         = Z$covariances,
        covariances.model   = covariances.model
    ), class = "cec");    
}

cec.interactive <- function(
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
        i <- 0    
        if (!is.matrix(centers)) centers <- init.centers(x, centers, centers.init)
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
            if (i == 0)
                desc = "(position of center means before first iteration)"      
            
            cat("Iterations:", Z$iterations, desc, "cost function:", Z$cost," \n ")
            
            plot(Z, ellipses = TRUE)
            
            if (readline) 
            {
                line = readline(prompt="Press <Enter> OR write number OR write 'q':");     
                lineint = suppressWarnings(as.integer(line))
                
                if (!is.na(lineint)) 
                {
                    i = i + lineint - 1
                    if (i < 0) 
                        i = -1        
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
            ignore = readline(prompt="Press <Enter>:")
        
        Z
    }
}