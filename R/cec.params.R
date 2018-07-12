# maps clustering type to int
resolve.type <- function(type)
{
    types <- c("covariance", "fixedr", "spherical", "diagonal", "eigenvalues", "mean", "all")
    match.arg(type, types)
}

# prepares clustering parameters for C function
create.cec.params.for.models <- function(k, n, type.arg, param.arg)
{
    models <- replicate(k, list())
    types <- vapply(type.arg, resolve.type, "")

    params <- NULL
    if (hasArg(param.arg))
        params <- param.arg
    if (length(types) == 1) {
        types <- rep(types, k)
        if (hasArg(param.arg)) {
            params <- rep(list(unlist(param.arg)), k)
            params <- params[!params %in% list(NULL, NA)]
        }
    }

    if (k != length(types))
        stop("Illegal argument: illegal length of \"type\" vector.")

    idx <- 0
    
    for (i in 1:length(types))
    {
        type = types[i]
        models[[i]]$type = type
        models[[i]]$params = list()
        if (type == resolve.type("covariance"))
        {
            idx <- idx + 1

            if (length(params) < idx)
                stop("Illegal argument: illegal param length.")

            cov <- params[[idx]]
            
            if (!is.array(cov)) stop("Illegal argument: illegal parameter for \"covariance\" type.")    
            if (ncol(cov) != n) stop("Illegal argument: illegal parameter for \"covariance\" type.")    
            if (nrow(cov) != n) stop("Illegal argument: illegal parameter for \"covariance\" type.")
            
            if (!try.chol(cov)) 
                stop("Illegal argument: illegal parameter for \"covariance\" type - matrix must be positive-definite.")
            
            cov.inv = solve(cov)
            models[[i]]$params <- list(cov = cov, cov.inv = cov.inv)
        }
        else if (type == resolve.type("fixed"))
        {
            idx <- idx + 1
            
            if (length(params) < idx)
                stop("Illegal argument: illegal param length.")
            
            r = params[[idx]]
            if (length(r) != 1) stop("Illegal argument: illegal parameter for \"fixedr\" type.")
            if (!is.numeric(r)) stop("Illegal argument: illegal parameter for \"fixedr\" type.")
            if (!r > 0)  stop("Illegal argument: illegal parameter for \"fixedr\" type.")
            models[[i]]$params = list(r = r)
        } 
        else if (type == resolve.type("eigenvalues"))
        {
            idx <- idx + 1
            
            if (length(params) < idx)
                stop("Illegal argument: illegal param length.")         
            
            evals <- params[[idx]]
            
            if (length(evals) != n) stop("Illegal argument: illegal parameter for \"eigenvalues\" type: invalid length.")
            if (!all(evals != 0)) stop("Illegal argument: illegal parameter for \"eigenvalues\" type: all values must be greater than 0.")
            models[[i]]$params = list(eigenvalues = sort(evals))
        }
        else if (type == resolve.type("mean"))
        {
            idx <- idx + 1
            
            if (length(params) < idx)
                stop("Illegal argument: illegal param length.")         
            
            mean <- params[[idx]]
            
            if (length(mean) != n) stop("Illegal argument: illegal parameter for \"mean\" type: invalid length.")
            models[[i]]$params = list(mean = mean)
        }
    }
    models
}

try.chol <- function(mat)
{
    ifelse("try-error" %in% class(try(chol(mat), silent=TRUE)), F, T)
}
