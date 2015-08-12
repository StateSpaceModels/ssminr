#'Build model in SSM
#'
#'This function build and compile a model in SSM. 
#' @param model_path character, full path to model directory.
#' @param pop_name character, name of the population.
#' @param data data.frame, must have a column named \code{date} that contains 
#' the observation dates (in \code{YYYY-MM-DD} format); and one or more numeric column(s) 
#' that contains the observed states and which are named according to the convention 
#' \code{x_obs}, where \code{x} is the name of the oberved state in the model. Missing
#' data are represented by \code{NA}.
#' @param start_date character, starting date of model integration (in \code{YYYY-MM-DD} format).
#' @param inputs list of inputs.
#' @param reactions list of reactions.
#' @param observations list of observations.
#' @inheritParams r2ssm
#' @return \code{ssm} object
#' @export
#' @aliases ssm
#' @import dplyr rjson
#' @importFrom plyr l_ply llply
#' @example inst/examples/SIR-example.r
new_ssm <- function(model_path, pop_name, data, start_date, inputs, reactions, observations) {

	# list directories
	if(!file.exists(model_path)){
		dir.create(model_path, recursive=TRUE)
	}

	start_date <- as.Date(start_date)

	wd <- setwd(model_path)
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

	# keep what you need
	data <- data %>% mutate(date=as.Date(date)) %>% filter(date > start_date)

	# write data
	data %>% mutate(date=as.character(date)) %>% write.csv(data_path,row.names=FALSE)

	time_series <- setdiff(names(data),"date")
	# link to ssm.json
	ssm_data <- plyr::llply(time_series, function(x) {
		list(name=x,require=list(path=data_path,fields=c("date",x)))		
	})
	
	# WRITE PRIORS ---------------------------------------------------------------------	
	# and return prior list for R
	priors <- plyr::llply(inputs, function(input){

		if(!is.null(input$prior)){

			input$prior %>% r2ssm_prior %>% rjson::toJSON(.) %>% write(file=file.path(dir_priors,paste0(input$name,".json")))

			if(input$prior$dist!="dirac"){
				prior <- input$prior
				prior$name <- input$name
				return(prior)	
			}
			
		}
	}) %>% remove_null

	# CREATE INPUTS ---------------------------------------------------------------------
	remainder_state <- NULL
	pop_size_theta <- NULL

	ssm_inputs <- plyr::llply(inputs,function(input) {

		if(!is.null(input$prior)){
			input$prior <- NULL
			input$require <- list(name=input$name,path=file.path(dir_priors,paste0(input$name,".json")))
		}

		# check tag
		if(input$tag=="remainder"){
			remainder_state <<- input$name
		}

		if(input$tag=="pop_size"){
			pop_size_theta <<- input$name
		}

		# remove value and tag
		input$value <- NULL
		input$tag <- NULL

		# remove all NULL
		input %>% remove_null %>% return
	})

	# remove remainder
	if(!is.null(remainder_state)){
		i_remainder <- which(get_name(ssm_inputs)!=remainder_state)
		ssm_inputs <- ssm_inputs[i_remainder]
	}

	# CREATE PROCESS ---------------------------------------------------------------------

	# unlist reaction if needed
	need_unlist <- sapply(reactions, function(x) x %>% names %>% is.null) 

	if(any(need_unlist)){

		index_unlist <- which(need_unlist)
		reactions_unlisted <- reactions[index_unlist] %>% unlist(recursive=FALSE, use.names=FALSE)
		reactions <- c(reactions[-index_unlist], reactions_unlisted)
		
	}


	# populations
	state_variables <- plyr::llply(reactions, function(x) {c(x$from,x$to)}) %>% unlist %>% unique
	ssm_populations <- r2ssm_populations(pop_name=pop_name, state_variables=state_variables, remainder=remainder_state, pop_size=pop_size_theta) 

	# CREATE OBSERVATIONS ---------------------------------------------------------------------

	# SSM currently needs the same start time for all observations.
	ssm_observations <- plyr::llply(observations, function(obs) {obs$start=as.character(start_date); return(obs)})

	# CREATE SDE ON INPUTS ---------------------------------------------------------------------

	# extract inputs with sde
	sde <- sapply(inputs, function(input) {
		# add input name and return
		x <- input$sde
		if(!is.null(x)){
			x$name <- input$name			
		}
		return(x)
	}) %>% remove_null
	
	if((n_sde <- length(sde))){

		# drift
		drift <- llply(sde, function(x) {

			tmp <- list(name=x$name, f=0)

			if(x$transformation!="none"){

				tmp$transformation <- switch(x$transformation,
					"log"=sprintf("log(%s)",x$name)
					)

			}

			return(tmp)
		})

		# dispersion matrix; only diagonal
		if(n_sde>1){

			input_sde <- get_name(sde)

			dispersion <- matrix(0, nrow=n_sde, ncol=n_sde, dimnames=list(input_sde,input_sde))

			diag(dispersion) <- sapply(sde, function(x) {x$volatility})

			colnames(dispersion) <- NULL
			dispersion <- as.data.frame(t(dispersion))
			colnames(dispersion) <- NULL		


		} else {

			dispersion <- list(list(sde[[1]]$volatility))

		}

		ssm_sde <- list(drift=drift, dispersion=dispersion)
		# cat(toJSON(ssm_sde))

	}

	# RESOURCES ---------------------------------------------------------------------

	# check which values are defined in input
	input_values <- sapply(inputs, function(input) {input$value})
	names(input_values) <- get_name(inputs)
	
	# extract theta from inputs: only inputs with a prior
	init_theta <- one_theta_sample_prior(priors) 

	# check if value is provided, if so set init
	theta_values <- input_values[names(init_theta)] %>% unlist
	init_theta[names(theta_values)] <- theta_values

	# default covmat
	init_covmat <- diag(init_theta/10)
	colnames(init_covmat) <- rownames(init_covmat) <- names(init_theta)

	ssm_theta <- r2ssm_resources(init_theta, init_covmat) 
	write(rjson::toJSON(ssm_theta),file=file.path(model_path,"theta.json"))

	# CREATE SSM FILES ---------------------------------------------------------------------

	ssm_json <- list(data=ssm_data, inputs=ssm_inputs, populations=ssm_populations, reactions=reactions, observations=ssm_observations) 
	if(n_sde){
		ssm_json$sde <- ssm_sde
	} 

	write(rjson::toJSON(ssm_json),file=file.path(model_path,"ssm.json"))

	# COMPILE MODEL ---------------------------------------------------------------------

	cmd <- sprintf("ssm -s %s/ssm.json",model_path)
	system(cmd)

	# RETURN SSM ---------------------------------------------------------------------
	
	setwd(wd)

	return(structure(list(
		model_path = model_path,
		pop_name = pop_name,
		state_variables = state_variables,
		theta = init_theta,
		covmat = init_covmat,
		summary = NULL,
		priors = priors,
		data = data,
		start_date = start_date,
		inputs = inputs,
		reactions = reactions,
		observations = observations),
	class="ssm"))
}


