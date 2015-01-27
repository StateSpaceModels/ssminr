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

#' @param theta a vector
#' @param covmat a matrix 
#' @name r2ssm
#' @aliases r2ssm_resources
r2ssm_resources <- function(theta, covmat){

	covmat_list <- as.list(as.data.frame(covmat))

	for(i in names(theta)){
		x <- covmat_list[[i]]
		x <- as.list(x)
		names(x) <- rownames(covmat)	
		covmat_list[[i]] <- x[x!=0] 
	}

	resources_formated <- list(
		list(name="values",description="initial values for the parameters",data=as.list(theta)),
		list(name="covariance",description="covariance matrix",data=covmat_list)
		)

	return(list(resources=resources_formated))

}




#' @param ssm_theta a list (usually a theta.json of SSM parsed by \code{\link[rjson]{fromJSON}})
#' @name r2ssm
#' @aliases ssm2r_resources
ssm2r_resources <- function(ssm_theta){

	resources <- ssm_theta$resources

	# theta
	theta <- resources[[1]]$data %>% unlist

	# covmat
	covmat <- resources[[2]]$data %>% plyr::ldply(as.data.frame, .id="row_names")
	rownames(covmat) <- covmat$row_names
	covmat$row_names <- NULL
	covmat[is.na(covmat)] <- 0


	if(length(resources)==3){
		# summary
		summary <- resources[[3]]$data %>% unlist
	} else {
		summary <- NULL
	}

	return(list(theta=theta, covmat=covmat, summary=summary))

}


#' @param pop_name character
#' @param state_variables character vector, name of all state variables 
#' @param remainder character, name of the state variable considered as a remainder.
#' @param pop_size character, name of the parameter corresponding to the population size. Required if \code{remainder} is provided.
#' @name r2ssm
#' @aliases r2ssm_populations
r2ssm_populations <- function(pop_name, state_variables, remainder=NULL, pop_size=NULL) {

	pop_formated <- list(name=pop_name, composition=state_variables)

	if(!is.null(remainder)){
		pop_formated$remainder <- list(name=remainder,pop_size=pop_size)
	}

	return(list(pop_formated))
}



r2ssm_args <- function(...) {

	# remove FALSE, NULL and replace TRUE by empty character
	args <- as.list(...) %>% clean_args

	# rename some arguments for SSM ('next' is a protected name in R)
	# Note: we don't use next, as it breaks package convention
	# that's ok as no input/output is erased by default
	# if("suffix"%in%names(args)){
	# 	names(args)[names(args)=="suffix"] <- "next"
	# }

	return(paste0("--",names(args)," ",args, collapse=" "))

}





