#'Convert R object into SSM object
#'
#'Internal functions to convert R object into SSM 
#' @param prior a list 
#' @importFrom plyr rename llply 
#' @name r2ssm
#' @aliases r2ssm_prior
r2ssm_prior <- function(prior) {

	if(prior$dist=="dirac"){
		return(list(name="dirac", distributionParameter=list(list(value=prior$args$value))))
	}

	dist_name <- switch(prior$dist,
		unif="uniform",
		truncnorm="normal"
		)

	dist_param <- plyr::rename(prior$args, c(a="lower",b="upper",min="lower",max="upper"),warn_missing=FALSE)
	dist_param[!is.finite(unlist(dist_param))] <- NULL

	dist_param <- plyr::llply(names(dist_param),function(param_name) {list(name=param_name,value=dist_param[[param_name]])})

	if( !is.null(prior$unit)){
		dist_param <- plyr::llply(dist_param,function(x) {c(x,unitCode=prior$unit)})
	}

	prior_formated <- list(name=dist_name, distributionParameter=dist_param)
	
	return(prior_formated)
}


