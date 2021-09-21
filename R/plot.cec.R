plot.cec <- function(x, col, cex = 0.5, pch = 19, cex.centers = 1, pch.centers = 8, ellipses.lwd = 4, ellipses = TRUE, model = TRUE, xlab, ylab, ...)
{
    if (ncol (x $ data) != 2 )
        stop("plotting available only for 2-dimensional data")
    
    if (!hasArg(col)) col = x$cluster
    
    if (!is.null(colnames(x$data))) 
    {
        xl <- colnames(x$data)[1]
        yl <- colnames(x$data)[2]
    } 
    else 
    {
        xl <- "x"
        yl <- "y"
    }
    
    if (hasArg(xlab)) xl <- xlab
    if (hasArg(ylab)) yl <- ylab
    
    plot(x$data, col=col, cex = cex, pch = pch, xlab = xl, ylab = yl, ...)
    
    if (model == T) {
        covs <- x$covariances.model
        means <- x$means.model
    } else {
        covs <- x$covariances
        means <- x$centers
    }
    
    points(x$means.model, cex = cex.centers, pch = pch.centers)   
    if (ellipses)
    {    
        for (i in 1:nrow(means))     
            if (! is.na(means[i, 1]))
            {         
                err = FALSE        
                tryCatch(
                    {
                        cov <- covs[[i]]
                        pts <- ellipse(means[i, ], cov)
                        lines(pts, lwd = ellipses.lwd)
                    },
                    finally = {})     
            }
    }  
}
