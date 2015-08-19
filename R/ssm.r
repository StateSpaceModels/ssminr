#'Build model in SSM
#'
#'This function build and compile a model in SSM. 
#' @param model_path character, full path to model directory.
#' @param pop character, name of the population.
#' @param data data.frame, must have a column named \code{date} that contains 
#' the observation dates (in \code{YYYY-MM-DD} format); and one or more numeric column(s) 
#' that contains the observed states and which are named according to the convention 
#' \code{x_obs}, where \code{x} is the name of the oberved state in the model. Missing
#' data are represented by \code{NA}.
#' @param start_date character, starting date of model integration (in \code{YYYY-MM-DD} format).
#' @param inputs list of inputs.
#' @param reactions list of reactions.
#' @param observations list of observations.
#' @param erlang_shapes named vector with state variables as names and shapes of the corresponding erlang distribution as values (optional).
#' @inheritParams r2ssm
#' @return \code{ssm} object
#' @export
#' @aliases ssm
#' @import dplyr rjson
#' @importFrom plyr l_ply llply dlply
#' @importFrom magrittr not
#' @example inst/examples/SIR-example.r
new_ssm <- function(model_path, pop, data, start_date, inputs, reactions, observations, erlang_shapes = NULL, states_in_SF = FALSE) {

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

	# check erlang_shapes
	erlang_shapes <- erlang_shapes[erlang_shapes > 1]

	if(length(erlang_shapes)==0){
		
		erlang_shapes <- NULL

	} else {

		# check that all erlang_shapes correspond to a state_variable
		state_variables <- get_state_variables(reactions)

		if(length(x <- setdiff(names(erlang_shapes), state_variables))){

			stop("The following elements of ", sQuote("erlang_shapes")," are not state variables: ", sQuote(x), call. = FALSE)
			
		}
	}


	# CREATE DATA ---------------------------------------------------------------------


	if(!all(c("date","time_series","value")%in%names(data))){
		stop("Missing column in ",sQuote("data")," see help(ssm)")
	}

	# keep what you need
	data <- data %>% mutate(date=as.Date(date), time_series = as.character(time_series)) %>% filter(date > start_date) 

	# write data
	ssm_data <- plyr::dlply(data, "time_series", function(df) {

		time_series <- first(df$time_series)
		data_path <- file.path(dir_data, paste0("ts_",time_series,".csv"))

		df %>% mutate(date=as.character(date)) %>% spread(time_series, value) %>% write.csv(data_path,row.names=FALSE)

		return(list(name=time_series,require=list(path=data_path,fields=c("date",time_series))))

	}) %>% unname
	
	# WRITE PRIORS ---------------------------------------------------------------------	
	# and return prior list for R

	priors <- plyr::llply(inputs, function(input){

		# if(0){
		# 	input <- find_element(inputs, "name", "N__pop_hcw")[[1]]
		# }

		if(!is.null(input$prior)){

			# print(input$name)

			input$prior %>% r2ssm_prior %>% rjson::toJSON(.) %>% write(file=file.path(dir_priors,paste0(input$name,".json")))

			if(input$prior$dist!="dirac"){
				prior <- input$prior
				prior$name <- input$name
				return(prior)	
			}
			
		}
	}) %>% remove_null

	# CREATE INPUTS ---------------------------------------------------------------------
	
	if(!is.null(erlang_shapes)){
		# make erlang
		inputs <- make_erlang_inputs(inputs, erlang_shapes) 
	}

	remainder_state <- find_element(inputs, "tag", "remainder") %>% get_name
	pop_size_theta <- find_element(inputs, "tag", "pop_size") %>% get_name

	ssm_inputs <- plyr::llply(inputs,function(input) {

		if(!is.null(input$prior)){
			input$prior <- NULL
			input$require <- list(name=input$name,path=file.path(dir_priors,paste0(input$name,".json")))
		}

		# remove value and tag
		input$value <- NULL
		input$tag <- NULL

		# remove all NULL
		input %>% remove_null %>% return
	})

	# remove remainder
	if(length(remainder_state)){
		i_remainder <- which(get_name(ssm_inputs)%in%remainder_state)
		ssm_inputs <- ssm_inputs[-i_remainder]
	}


	# CREATE REACTIONS ---------------------------------------------------------------------

	# unlist reaction if needed
	need_split <- sapply(reactions, function(x) length(x$to)>1) 

	if(any(need_split)){

		index_split <- which(need_split)
		for(i in index_split){

			r <- reactions[[i]]
			reactions <- c(reactions,reaction_split(from=r$from, split_to=r$to, description=r$description, rate=r$rate, accumulators=r$accumulators, keywords=r$keywords))

		}

		reactions <- reactions[-index_split]
	}

	if(!is.null(erlang_shapes)){
		# make erlang
		reactions <- make_erlang_reactions(reactions, erlang_shapes) 
	}

	if(any(need_split)){

		i_split <- get_element(reactions, "split") %>% sapply(is.null) %>% magrittr::not() %>% which
		
		for(i in i_split){

			reactions[[i]]$rate <- sprintf("(%s)*(%s)", reactions[[i]]$split, reactions[[i]]$rate)
			reactions[[i]]$split <- NULL

		}

	}

	# divide the rate of linear reactions to compensate for density-dependence in SSM
	# TODO: implement a special function or a way to tell SSM how to deal with linear reactions
	i_linear <- get_element(reactions, "keywords") %>% sapply(function(x) "linear"%in%x) %>% which
	
	for(i in i_linear){
		reactions[[i]]$rate <- sprintf("(%s)/(%s)", reactions[[i]]$rate, reactions[[i]]$from)		
	}


	# ensure positivity of from compartment using a special function (1_{from > 0})
	# using states in special function break compilation as special functions aren't derived 
	# during the diffusion approximation => flag and use dirty fix for compilation (see at the end)
	# NOTE: this fix prevent the use of sde approx and kalman medthods in SSM
	# TODO: SSM should be able to deal with state within special functions	

	i_while <- get_element(reactions, "keywords") %>% sapply(function(x) "while_from_is_positive"%in%x) %>% which
	
	if(length(i_while)){
		# flag to use hack form compilation
		states_in_SF <- TRUE
	}

	for(i in i_while){
		reactions[[i]]$rate <- sprintf("(%s)*heaviside(%s - 1)", reactions[[i]]$rate, reactions[[i]]$from)		
	}


	# CREATE POPULATIONS ---------------------------------------------------------------------

	state_variables <- get_state_variables(reactions)
	ssm_populations <- r2ssm_populations(pop=pop, state_variables=state_variables, remainder=remainder_state, pop_size=pop_size_theta) 

	## replace remainder in reaction rates if it is present
	# TODO: SSM should deal with that automatically
	if(length(remainder_state)){

		replace_remainder <- sapply(ssm_populations, function(ssm_population) {
			sprintf("(%s - %s)", ssm_population$remainder$pop_size, ssm_population$composition %>% setdiff(ssm_population$remainder$name) %>% paste(collapse=" - "))

		})

		names(replace_remainder) <- remainder_state

		reactions <- llply(reactions, function(reaction) {

			reaction$rate <- reaction$rate %>% str_replace_all(replace_remainder)

			return(reaction)
		})


	}

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

	if(states_in_SF){
	# TEMP HACK FOR FIXING ISSUE WITH SPECIAL FUNCTION THAT NEED TO BE DERIVED WHEN THEY DEPEND ON STATE VARIABLES
    # the problem is with jac.c, the hack consists in using the jac.c generated by the same model without state_variable dependencies in special functions
    # this file will be created in ../bin/C/templates/jac.c

	# 1- replace all reactions rates by 1
		ssm_json_hack <- ssm_json
		ssm_json_hack$reactions <- llply(ssm_json_hack$reactions, function(reaction) {

			reaction$rate <- "1"

			return(reaction)

		})

	# 2- compile silently the model without rates
		write(rjson::toJSON(ssm_json_hack),file=file.path(model_path,"ssm.json"))
		system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

	# 3- save jac.c in /tmp
		file.copy(file.path(model_path, "bin", "C", "templates", "jac.c"), tempdir(), overwrite = TRUE)

	# 4- compile silently the model with rates
		write(rjson::toJSON(ssm_json),file=file.path(model_path,"ssm.json"))
		system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
		
	# 5- replace jac.c by tmp/jac.c
		file.copy(file.path(tempdir(), "jac.c"), file.path(model_path, "bin", "C", "templates"), overwrite = TRUE)

	# 6- recompile manually the model
		cmd <- sprintf("cd %s/bin/C/templates; make clean; make; make install",model_path)
		system(cmd)

		message("Having states in special functions prevents the use of sde approx and kalman medthods in SSM")

	} else {

		system(cmd)

	}


	# RETURN SSM ---------------------------------------------------------------------

	setwd(wd)

	return(structure(list(
		model_path = model_path,
		pop = pop,
		state_variables = state_variables,
		theta = init_theta,
		covmat = init_covmat,
		summary = NULL,
		priors = priors,
		data = data,
		start_date = start_date,
		inputs = inputs,
		reactions = reactions,
		observations = observations,
		erlang_shapes=erlang_shapes),
	class="ssm"))
}


