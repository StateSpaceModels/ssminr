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

	return(list(name=x$obs, distribution="discretized_normal",mean=sprintf("%s * (%s)", reporting, x$state), sd = sprintf("sqrt(%s * ( 1.0 - %s ) * (%s) + pow(%s * %s * (%s),2))", reporting, reporting,  x$state, reporting, overdispersion, x$state)))

}