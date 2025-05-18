#' tsclean2
#`
#` Extends the forecast package's tsoutliers function to replace missing
#' values as well, using a user-specified function from the imputeTS package.
#' @param x The vector to be cleansed
#' @param replace.missing If TRUE, it not only replaces outliers, but also replaces missing values
#' @param iterate The number of iterations required; pass-through to tsoutliers.
#' @param lambda  Box-Cox transformation parameter. If lambda="auto", then a transformation 
#' is automatically selected using BoxCox.lambda. The transformation is ignored if NULL. 
#' Otherwise, data transformed before model is estimated. Pass-through to tsoutliers.
#' ... Other parameters to be passed to method
#' @param method See imputeTS documentation for interpolation methods
#' @return  A vector with outliers replaced by imputed values

tsclean2 <- function (x, 
											replace.missing = TRUE, 
											iterate = 2, 
											lambda = NULL, 
											method = na_interpolation, ...) 
{
	require(forecast)
	require(imputeTS)
	
	outliers <- tsoutliers(x, iterate = iterate, lambda = lambda)
	x[outliers$index] <- outliers$replacements
	if (replace.missing) {
		x <- method(x, ...)
	}
	return(x)
}

#' clean
#'
#' Identify and replace outliers, including NAs and long runs of zero values,
#' in an integer-valued time series
#' @param x The vector to be cleansed
#' @param cutoff  The cutoff p-value used to identify "too long" a run of zero values
#' @param randomize   Randomize the replacement values (using a Poisson dist'n)
#' @param method  The  method used to calculate replacement values
#' @param nonnegative_only TRUE if negative return values replaced by zeros
#' @return A vector with values replaced
#' @examples 
#' x <- rpois(50,5)
#' x[31:35] <- 0
#' z <- clean(x, method = na_kalman, randomize=TRUE)
#' print(z[31:35])
#' @export

clean <- function(x, 
									cutoff = 0.999, 
									randomize = FALSE, 
									method = na_interpolation,
									nonnegative_only = TRUE) {
	
	require(tseries)
	
	y <- x
	y[y>0] <- 1
	run_info <- rle(y)
	runs_pvalue <- runs.test(as.factor(y))$p.value
	
	if (runs_pvalue < 1 - cutoff) {
		
		maxlr <- qgeom(cutoff, nones/(nzeros + nones)) - 1
		
		run_info$starts <- 1 + c(0, cumsum(run_info$lengths))
		for (i in seq_along(run_info$values)) {
			if (run_info$values[i] == 0 & run_info$lengths[i] > maxlr) {
				x[run_info$starts[i]:(run_info$starts[i] + run_info$lengths[i]-1)] <- NA
			}
		}
	}
	
	z <- (tsclean2(sqrt(x), replace.missing = TRUE, method = method))^2
	
	if (randomize) {
		outliers <- which(round(z) != x | is.na(x))
		if (length(outliers) > 0) {
			z[outliers] <- rpois(length(outliers), z[outliers])
		}
	}
	
	if (nonnegative_only) {
		return(pmax(0,round(z)))
	} 
	return(round(z))
}
