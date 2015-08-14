pop_name <- function(x, pop) {

	return(sprintf("%s__pop_%s" ,x, pop))

}

define_NGM <- function(S_pop, I_pop, prior, value, infectious_states, infectious_period, par_pop_size, base_name_NGM = "R", base_name_WAIFW = "beta") {

	if(0){
		S_pop <- pop
		I_pop <- pop
		prior <- unif(0,50)
		value <- 1
		infectious_states <- c("I_R", "I_D")
		infectious_period <- "d_infectious"
		par_pop_size <- "N"
		base_name_NGM <- "R"
	}

	#
	WAIFW <- expand.grid(S=S_pop, I=I_pop) %>% mutate(R=sprintf("%s_S%s_I%s",base_name_NGM,S,I), beta=sprintf("%s_S%s_I%s",base_name_WAIFW,S,I), infectious_period = pop_name(infectious_period,I), S_pop_size = pop_name(par_pop_size, S))

	# define inputs
	NGM_inputs <- dlply(WAIFW, c("R"), function(df) {

		list(
			input(name=df$R, description="basic reproduction number", prior=prior, value= value),
			input(name=df$beta, description="effective contact rate", transformation=sprintf("(%s)/(%s)", df$R, df$infectious_period))
			)

	}) %>% unlist(recursive=FALSE) %>% unname
	

	# define force of infection
	force_of_infection <- dlply(WAIFW, "S", function(df) {

		infectious_in_each_pop <- sapply(df$I, function(x) pop_name(infectious_states, x) %>% paste(collapse=" + ")) 

		sprintf("%s * (%s)",df$beta, infectious_in_each_pop) %>% paste(collapse=" + ") %>% sprintf("(%s) / (%s)", ., first(df$S_pop_size))

	}) %>% unlist(recursive=FALSE)

	return(list(inputs=NGM_inputs, force_of_infection=force_of_infection))

}

add_pop_to_input <- function(input, pop, names_var_pop) {


	# if(input$name == "N"){
	# 	browser()
	# }

	for(name_var_pop in names_var_pop){

		input$name <- str_replace_all(input$name, sprintf("\\b%s\\b",name_var_pop), pop_name(name_var_pop, pop))							

		if(!is.null(input$transformation)){
			input$transformation <- str_replace_all(input$transformation, sprintf("\\b%s\\b",name_var_pop), pop_name(name_var_pop, pop))							
		}

	}	

	# pick prior if necessary
	if(!is.null(names(input$prior)) && pop%in%names(input$prior)){
		input$prior <- input$prior[[pop]]
	}

	# pick value if necessary
	if(!is.null(names(input$value)) && pop%in%names(input$value)){
		input$value <- input$value[[pop]]
	}

	if(length(input$value)>1){
		stop("Too many input value for ", sQuote(input$name))
	}

	
	return(input)

}


add_pop_to_reaction <- function(reaction, pop, names_var_pop) {

	for(names_pop_input in names_var_pop){

		# change from
		reaction$from <- str_replace_all(reaction$from, sprintf("\\b%s\\b",names_pop_input), pop_name(names_pop_input, pop))				
		
		# change to
		names_to <- names(reaction$to)
		reaction$to <- str_replace_all(reaction$to, sprintf("\\b%s\\b",names_pop_input), pop_name(names_pop_input, pop))				
		if(!is.null(names_to)){
			# split reaction
			names(reaction$to) <- str_replace_all(names_to, sprintf("\\b%s\\b",names_pop_input), pop_name(names_pop_input, pop))				

		}
		
		# change rate
		if(!is.null(names(reaction$rate))){
			reaction$rate <- reaction$rate[pop]
		}

		if(length(reaction$rate)!=1){
			stop("Wrong specification of reaction rate: ", sQuote(reaction$rate))
		}

		reaction$rate <- str_replace_all(reaction$rate, sprintf("\\b%s\\b",names_pop_input), pop_name(names_pop_input, pop))				
		

		if(!is.null(reaction$accumulators)){
			# change accumulators
			reaction$accumulators[[1]] <- str_replace_all(reaction$accumulators[[1]], sprintf("\\b%s\\b",names_pop_input), pop_name(names_pop_input, pop))				
		}

		if(!is.null(reaction$white_noise)){
			# change white noise
			reaction$white_noise$name <- str_replace_all(reaction$white_noise$name, sprintf("\\b%s\\b",names_pop_input), pop_name(names_pop_input, pop))				
			reaction$white_noise$sd <- str_replace_all(reaction$white_noise$sd, sprintf("\\b%s\\b",names_pop_input), pop_name(names_pop_input, pop))				
		}

	}	

	return(reaction)

}




add_pop_to_observation <- function(observation, pop, names_var_pop) {

	for(name_var_pop in names_var_pop){

		observation$name <- str_replace_all(observation$name, sprintf("\\b%s\\b",name_var_pop), pop_name(name_var_pop, pop))							
		observation$mean <- str_replace_all(observation$mean, sprintf("\\b%s\\b",name_var_pop), pop_name(name_var_pop, pop))							
		
		if(!is.null(observation$sd)){
			observation$sd <- str_replace_all(observation$sd, sprintf("\\b%s\\b",name_var_pop), pop_name(name_var_pop, pop))										
		}

	}	

	return(observation)

}




make_pop_struct <- function(pop, inputs, reactions, observations, names_shared_inputs=NULL, erlang_shapes=NULL) {

	if(0){

		pop <- pop
		inputs <- SEIRD_inputs
		reactions <- SEIRD_reactions
		observations <- SEIRD_observations
		names_shared_inputs <- names_shared_inputs
		erlang_shapes <- Erlang_shapes

	}

	# struct input
	names_var_pop <- setdiff(get_name(inputs), names_shared_inputs)

	inputs_pop <- llply(inputs, function(input) {

		input_pop <- llply(pop, add_pop_to_input, input=input, names_var_pop=names_var_pop)

	}) %>% unlist(recursive = FALSE) %>% unique

	# inputs_pop %>% get_name %>% sort

	# struct reactions
	## add accumulators to names_var_pop
	names_accumulators_pop <- get_element(reactions, "accumulators") %>% unlist %>% unique %>% setdiff(names_shared_inputs)
	names_var_pop <- c(names_var_pop, names_accumulators_pop)

	reactions_pop <- llply(reactions, function(reaction) {

		reaction_pop <- llply(pop, add_pop_to_reaction, reaction=reaction, names_var_pop=names_var_pop)

	}) %>% unlist(recursive = FALSE) %>% unique

	# struct observation
	## add obs names to names_var_pop
	names_obs_pop <- get_name(observations) %>% unlist %>% unique %>% setdiff(names_shared_inputs)
	names_var_pop <- c(names_var_pop, names_obs_pop)

	observations_pop <- llply(observations, function(observation) {

		observation_pop <- llply(pop, add_pop_to_observation, observation=observation, names_var_pop=names_var_pop)

	}) %>% unlist(recursive = FALSE) %>% unique


	# struct erlang
	if(is.list(erlang_shapes) && is.null(names(erlang_shapes))){
		stop("List of erlang_shapes need pop names", sQuote(erlang_shapes))
	} else if (is.atomic(erlang_shapes)){
		erlang_shapes <- rep(list(erlang_shapes), length(pop))
		names(erlang_shapes) <- pop
	}

	df_erlang_shapes <- ldply(erlang_shapes, function(x) {data_frame(state = names(x), shape = x)}, .id="pop")

	df_erlang_shapes_pop <- df_erlang_shapes %>% filter(state %in% names_var_pop) %>% mutate(state = pop_name(state, pop)) %>% select(-pop)
	df_erlang_shapes_shared <- df_erlang_shapes %>% filter(!state %in% names_var_pop) %>% select(-pop) %>% distinct
	df_erlang_shapes_pop <- bind_rows(df_erlang_shapes_pop, df_erlang_shapes_shared)

	erlang_shapes_pop <- df_erlang_shapes_pop$shape
	names(erlang_shapes_pop) <- df_erlang_shapes_pop$state

	return(list(inputs=inputs_pop, reactions=reactions_pop, observations=observations_pop, erlang_shapes=erlang_shapes_pop))

}