#'Build model in SSM
#'
#'This function build and compile the model in SSM. 
#' @param model character, full path to model directory.
#' @param pop_name character, name of the population.
#' @param data data.frame, must have a column named \code{date} that contains 
#' the observation dates (in \code{YYYY-MM-DD} format); and one or more numeric column(s) 
#' that contains the observed states and which are named according to the convention 
#' \code{x_obs}, where \code{x} is the name of the oberved state in the model. Missing
#' data are represented by \code{NA}.
#' @param start_date character, starting date of model integration (in \code{YYYY-MM-DD} format).
#' @param priors list of priors.
#' @param inputs list of inputs.
#' @param reactions list of reactions.
#' @param observations list of observations.
#' @param diffed character, parameter names that follow a diffusion
#' @inheritParams r2ssm
#' @export
#' @import dplyr rjson
#' @importFrom plyr l_ply llply
build_ssm <- function(model, pop_name, data, start_date, priors, inputs, reactions, observations, remainder=NULL, pop_size=NULL, diffed=NULL) {

	# list directories
	if(!file.exists(model)){
		dir.create(model, recursive=TRUE)
	}

	start_date <- as.Date(start_date)


	wd <- setwd(model)
	# relative directories
	dir_data <- "data"
	dir_priors <- "priors"

	# create directories
	dir_list <- c(dir_data,dir_priors)

	plyr::l_ply(dir_list,function(dir) {
		if(!file.exists(dir)){
			dir.create(dir)
		}
	})

	# CREATE DATA ---------------------------------------------------------------------

	data_path <- file.path(dir_data,paste0("data_",pop_name,".csv"))

	# 
	data <- data %>% mutate(date=as.Date(date)) %>% filter(date > start_date) %>% mutate(date=as.character(date))

	# write data
	write.csv(data,data_path,row.names=FALSE)

	time_series <- setdiff(names(data),"date")
	# link to ssm.json
	data_formated <- plyr::llply(time_series, function(x) {
		list(name=x,require=list(path=data_path,fields=c("date",x)))		
	})
	
	# CREATE PRIORS ---------------------------------------------------------------------

	# write prior
	for(theta_name in names(priors)){

		x <- priors[[theta_name]]
		x <- r2ssm_prior(x)
		write(rjson::toJSON(x),file=file.path(dir_priors,paste0(theta_name,".json")))

	}

	# CREATE INPUTS ---------------------------------------------------------------------

	inputs_formated <- plyr::llply(names(inputs),function(input_name) {

		input <- c(name=input_name, inputs[[input_name]])

		if(!is.null(input$prior)){
			prior_name <- input$prior
			input$prior <- NULL
			input$require <- list(name=prior_name,path=file.path(dir_priors,paste0(prior_name,".json")))
		}

		return(input)
	})

	# CREATE PROCESS ---------------------------------------------------------------------

	# populations
	state_variables <- plyr::llply(reactions, function(x) {c(x$from,x$to)}) %>% unlist %>% unique
	populations_formated <- r2ssm_populations(pop_name=pop_name, state_variables=state_variables, remainder=remainder, pop_size=pop_size) 

	# CREATE OBSERVATIONS ---------------------------------------------------------------------

	observations_formated <- plyr::llply(observations, r2ssm_observation, start_date=start_date)

	# CREATE diffed ---------------------------------------------------------------------
	diffed_theta <- diffed
	
	if(!is.null(diffed)){

		dispersion <- matrix(0, nrow=length(diffed_theta), ncol=length(diffed_theta), dimnames=list(diffed_theta,diffed_theta))

		diag(dispersion) <- "vol"
		
		drift <- llply(diffed_theta, function(x) {

			if(str_detect(x,"beta")){				
				return(list(name=x, f=0, transformation=sprintf("log(%s)",x)))
			}

		})

		colnames(dispersion) <- NULL
		dispersion <- as.data.frame(t(dispersion))
		colnames(dispersion) <- NULL		

		if(all(dim(dispersion)==c(1,1))){
			dispersion <- list(list("vol"))
		}

		diffed_formated <- list(drift=drift, dispersion=dispersion)
		# cat(toJSON(diffed_formated))

	}

	# RESOURCES ---------------------------------------------------------------------

	# read prior slot of ssm_inputs for theta_names
	# look at theta_priors, remove if dirac
	theta_dist <- unlist(llply(priors, function(x) {x$dist}))
	theta_estimated_names <- names(priors)[theta_dist!="dirac"]

	init_theta <- sample_from_prior(priors, theta_estimated_names) 
	names(init_theta) <- theta_estimated_names

	init_covmat <- diag(init_theta/10)
	colnames(init_covmat) <- rownames(init_covmat) <- names(init_theta)

	resources_formated <- r2ssm_resources(init_theta, init_covmat) 

	# CREATE SSM FILES ---------------------------------------------------------------------

	ssm_json <- list(data=data_formated, inputs=inputs_formated, populations=populations_formated, reactions=reactions, observations=observations_formated) 
	if(!is.null(diffed)){
		ssm_json$sde <- diffed_formated
	} 

	theta_json <- list(resources=resources_formated)

	write(rjson::toJSON(ssm_json),file=file.path(model,"ssm.json"))
	write(rjson::toJSON(theta_json),file=file.path(model,"theta.json"))

	# COMPILE MODEL ---------------------------------------------------------------------

	cmd <- sprintf("ssm -s %s/ssm.json",model)
	system(cmd)

	# SAVE MODEL IN R ---------------------------------------------------------------------
	
	ssm_model <- list(data=data, inputs=inputs, state_variables=state_variables, theta=list(init=init_theta, covmat=init_covmat), theta_priors=priors)
	saveRDS(ssm_model,file.path(model,"ssm_model.rds"))

	setwd(wd)


}



