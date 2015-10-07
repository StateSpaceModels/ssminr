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

	names_theta <- names(theta)
	add_theta <- setdiff(names_theta, rownames(covmat))

	if(length(add_theta)){
		# add missing theta to covmat
		# for instance when mcmc(ode) -> kmcmc(sde) and one wants to estimate volatility
		n_theta <- length(theta)
		new_covmat <- matrix(0, nrow=n_theta, ncol=n_theta, dimnames=list(names_theta,names_theta))
		new_covmat[rownames(covmat), colnames(covmat)] <- as.matrix(covmat)
		new_covmat[add_theta,add_theta] <- theta[add_theta]/10
		covmat <- new_covmat
	}

	covmat_list <- as.list(as.data.frame(covmat))

	for(i in names(covmat_list)){
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


r2ssm_one_population <- function(pop, state_variables, remainder=NULL, pop_size=NULL) {

	if(length(pop)!=1){
		stop("More than one pop")
	}

	pop_formated <- list(name=pop, composition=state_variables)

	if(!is.null(remainder)){

		if(length(remainder) > 1){
			stop("More than one remainder_state ", sQuote(remainder))	
		}

		if(length(pop_size) != 1){
			stop("Only one pop_size required. Pop_size found: ", sQuote(pop_size))	
		}


		pop_formated$remainder <- list(name=remainder,pop_size=pop_size)
	}

	return(pop_formated)
}


#' @param pop character
#' @param state_variables character vector, name of all state variables 
#' @param remainder character, name of the state variable considered as a remainder.
#' @param pop_size character, name of the parameter corresponding to the population size. Required if \code{remainder} is provided.
#' @name r2ssm
#' @aliases r2ssm_populations
r2ssm_populations <- function(pop, state_variables, remainder=NULL, pop_size=NULL) {

	if(length(pop)==1){

		pop_formated <- list(r2ssm_one_population(pop, state_variables, remainder, pop_size))

	} else {

		pop_formated <- llply(pop, function(pop_x) {

			state_variables_x <- state_variables[str_detect(state_variables, sprintf("pop_%s", pop_x))]

			if(!is.null(remainder))	{
				# search remainder of the pop
				remainder_x <- remainder[str_detect(remainder, sprintf("pop_%s", pop_x))]			
			} else {
				remainder_x <- NULL
			}

			if(!is.null(pop_size))	{
				# search pop_size of the pop
				pop_size_x <- pop_size[str_detect(pop_size, sprintf("pop_%s", pop_x))]			
			} else {
				pop_size_x <- NULL
			}

			r2ssm_one_population(pop_x, state_variables_x, remainder_x, pop_size_x) %>% return

		})
		
	}

	return(pop_formated)
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

	# TODO: better patch
	if("iter"%in%names(args)){
		args[["iter"]] <- args[["iter"]] %>% as.numeric %>% format(scientific=FALSE)
	}

	return(paste0("--",names(args)," ",args, collapse=" "))

}





