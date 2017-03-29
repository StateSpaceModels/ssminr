#' @import stringr
set_obs_variable <- function(obs=NULL, state=NULL) {

	if(is.null(obs) && is.null(state)){
		stop("At least one of obs or state must be provided")
	}

	if(is.null(obs) && !is.null(state)){
		# only state provided
		obs <- sprintf("%s_obs",state)
	} else if(is.null(state)){
		# only obs provided
		state <- str_replace(obs, "_obs", "")
	} 

	return(list(obs=obs, state=state))

}


#'Observation processes for SSM
#'
#'SSM currently accept the following observation processes.
#' @param obs character, name of the observation in the data
#' @param  state character, name of the observed state in the model
#' @param  reporting character, name of the reporting rate parameter
#' @export
#' @name obs
#' @aliases poisson_obs
poisson_obs <- function(obs=NULL, state=NULL, reporting) {

	x <- set_obs_variable(obs, state)

	return(list(name=x$obs, distribution="poisson",mean=sprintf("%s * (%s)", reporting, x$state)))

}

#' @param  overdispersion character, name of the overdispersion parameter
#' @export
#' @name obs
#' @aliases discretized_normal_obs
discretized_normal_obs <- function(obs=NULL, state=NULL, reporting, overdispersion) {

	x <- set_obs_variable(obs, state)

	return(list(name=x$obs, distribution="discretized_normal",mean=sprintf("%s * (%s)", reporting, x$state), sd = sprintf("sqrt(%s * (%s) + %s * pow(%s * (%s),2))", reporting,  x$state, overdispersion, reporting, x$state)))

}


#' @param  x_obs character, name of the column with the number of successes in the data
#' @param  n_obs character, name of the column with the sample size in the data
#' @param  p character, probability of success, usually a combination of the states and parameters of the model (e.g. \code{p="I/N"})
#' @export
#' @name obs
#' @aliases binomial_obs
binomial_obs <- function(x_obs=NULL, n_obs=NULL, p=NULL) {

	return(list(name=x_obs, distribution="binomial", n = n_obs, p = p))

}


