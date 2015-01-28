#'Some priors distributions
#'
#'Distributions supported by SSM to define priors
#' @inheritParams stats::dunif
#' @export
#' @seealso \code{\link{one_theta_sample_prior}}
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
truncnorm <- function(a=-Inf, b=Inf, mean, sd){

	return(list(dist="truncnorm",args=list(a = a, b = b, mean=mean, sd=sd)))

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
#' @param theta_names character, names of the parameters sampled. If \code{NULL} (default) all parameters are sampled.
#' @export
one_theta_sample_prior <- function(priors, theta_names=NULL) {

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


#'Sample prior
#'
#'Generate one sample \code{theta} of \code{ssm} from the prior distribution
#' @inheritParams call_ssm
#' @export
#' @seealso \code{\link{one_theta_sample_prior}}
#' @return a \code{ssm} object
one_ssm_sample_prior <- function(ssm) {

	if(!inherits(ssm,"ssm")){
		stop(sQuote("ssm"),"is not an object of class ssm")
	}

	ssm$theta <- one_theta_sample_prior(ssm$prior)

	return(ssm)
}


#'Sample prior
#'
#'Generate one or more \code{ssm} objects with \code{theta} sampled from the prior distribution.
#' @param n numeric, sample size
#' @param method character, method used to sample from prior distribution:
#' \itemize{
#' 	\item \code{"random"} randomly sample from the prior distribution
#' }
#' @inheritParams call_ssm
#' @export
#' @import dplyr
#' @seealso \code{\link{one_theta_sample_prior}} \code{\link{one_ssm_sample_prior}}
#' @return if \code{n==1}: a \code{ssm} object, otherwise a \code{\link{tbl}} object with two columns: \code{id} (sampel index) and \code{ssm} (corresponding \code{ssm} object)
sample_prior <- function(ssm, n, method=c("random")) {	

	stopifnot(n>0)

	method <- match.arg(method)

	sampler <- switch(method,
		"random"="one_ssm_sample_prior" 
		)

	if(n==1){
		do.call(sampler, list(ssm=ssm))
	} else {
		data_frame(id=1:n) %>% group_by(id) %>% do(ssm=do.call(sampler, list(ssm=ssm)))
	}
	
}





