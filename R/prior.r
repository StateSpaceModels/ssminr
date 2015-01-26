#'Some priors distributions
#'
#'Distributions supported by SSM to define priors
#' @inheritParams stats::dunif
#' @export
#' @seealso \code{\link{sample_from_prior}}
#' @name prior
#' @aliases unif
unif <- function(min, max){

	return(list(dist="unif",args=list(min = min, max = max)))

}

#' @inheritParams truncnorm::dtruncnorm
#' @export
#' @import truncnorm
#' @name prior
#' @aliases truncnorm
truncnorm <- function(a, b, mean, sd){

	return(list(dist="truncnorm",args=list(a = 0, b = 1, mean=0.6, sd=0.1)))

}

#' @param value numeric, dirac value
#' @export
#' @name prior
#' @aliases dirac
dirac <- function(value){

	return(list(dist="dirac",args=list(value=value)))

}


#'Sample from prior
#'
#'Generate a sample from prior density distribution.
#' @param priors a list of priors
#' @param  theta_names character, names of the parameters sampled. If \code{NULL} (default) all parameters are sampled.
#' @export
sample_from_prior <- function(priors, theta_names=NULL) {

	names(priors) <- get_name(priors)

	# restrict priors to theta_names if provided
	if(!is.null(theta_names)){
		priors <- priors[intersect(names(priors), theta_names)]
	}


	theta_sample <- sapply(priors, function(prior) {

		rprior <- paste0("r",prior$dist)
		prior_args <- c(prior$args,n=1)
		return(do.call(rprior,prior_args))

	})	

	return(theta_sample)
}
