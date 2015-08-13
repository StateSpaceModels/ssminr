make_erlang_reactions <- function(reactions, erlang_shapes) {

	erlang_states <- names(erlang_shapes)

	for(i in seq_along(erlang_shapes)){

		erlang_state <- erlang_states[i]
		erlang_shape <- erlang_shapes[[i]]

    	# change FROM
		i_reaction <- sapply(reactions, function(r) {r$from==erlang_state}) %>% which

		for(i_reaction_from in i_reaction){

			reaction <- reactions[i_reaction_from]

	        # if reaction is a "transmission" then just duplicate reaction and change "from"
			if(!is.null(reaction[[1]]$keywords) && c("transmission", "linear")%in%unlist(reaction[[1]]$keywords)){

				reaction_erlang <- rep(reaction, erlang_shape)

				reaction_erlang <- llply(seq_along(reaction_erlang), function(i) {

					reaction <- reaction_erlang[[i]]
					reaction$from <- paste0(erlang_state, i)
					return(reaction)
				})

			} else {

       			 # erlang waiting time
				reaction_rate <- reaction[[1]]$rate

       			 # mutiply rate by shape
				if(str_detect(reaction_rate, "correct_rate")){
					reaction_rate <- str_replace(reaction_rate, fixed("correct_rate("), paste0("correct_rate(",erlang_shape,"*"))

				} else {
					reaction_rate <- sprintf("%s*(%s)",erlang_shape,reaction_rate)
				}

				reaction[[1]]$rate <- reaction_rate

				reaction_erlang <- rep(reaction, erlang_shape)


				if(!is.null(reaction[[1]]$accumulators) & erlang_shape > 1){

          			# only last transition has an accumulator
					for(k in seq_len(erlang_shape-1)){
						reaction_erlang[[k]]$accumulators <- NULL
					}

				}

				if(!is.null(reaction[[1]]$split) & erlang_shape > 1){

          			# only last transition has a split
					for(k in seq_len(erlang_shape-1)){
						reaction_erlang[[k]]$split <- NULL
					}

				}

				reaction_erlang <- llply(seq_along(reaction_erlang), function(i) {

					reaction <- reaction_erlang[[i]]
					reaction$from <- paste0(erlang_state, i)
					reaction$to <- paste0(erlang_state, i+1)
					return(reaction)
				})

				reaction_erlang[[length(reaction_erlang)]]$to <- reaction[[1]]$to


				if("split"%in%names(reaction[[1]]) && !is.null(reaction[[1]]$keywords) && !"split1"%in%(reaction[[1]]$keywords %>% unlist)){

					# keep only last
					reaction_erlang <- reaction_erlang[length(reaction_erlang)]

				}

			}

			reactions <- c(reactions,reaction_erlang)

		}

		reactions[i_reaction] <- NULL
		
   		# change TO
		i_reaction <- sapply(reactions, function(r) {r$to==erlang_state}) %>% which
		for(i_reaction_to in i_reaction){
			reactions[[i_reaction_to]]$to <- paste0(erlang_state,1)
		}

		# change rate
		i_reaction <- sapply(reactions, function(r) {str_detect(r$rate,erlang_state)}) %>% which
		for(i_reaction_rate in i_reaction){
			reactions[[i_reaction_rate]]$rate <- str_replace_all(reactions[[i_reaction_rate]]$rate, sprintf("\\b%s\\b",erlang_state), paste0(erlang_state,1:erlang_shape, collapse=" + ") %>% protect)
		}


	}

	return(reactions)

}



make_erlang_inputs <- function(inputs, erlang_shapes) {

	# restrict to all erlang_states that are in inputs
	erlang_shapes <- erlang_shapes[intersect(names(erlang_shapes), get_name(inputs))]	

	erlang_states <- names(erlang_shapes)

	names(inputs) <- get_name(inputs)

	for(erlang_state in erlang_states){

		erlang_shape <- erlang_shapes[erlang_state]
		input <- inputs[erlang_state]

		if(!is.null(input[[1]]$prior)){

			# keep origianl input for prior and generate erlang input with transformation
			new_input <- input
			
			new_input[[1]]$transformation <- sprintf("(%s)/(%s)", input[[1]]$name, erlang_shape)
			new_input[[1]]$prior <- NULL
			
			# rescale initial value
			if(!is.null(input[[1]]$value)){
				new_input[[1]]$value <- sprintf("(%s)/(%s)", input[[1]]$value, erlang_shape)				
			}

			input_erlang <- rep(new_input, erlang_shape)

		} else {

			input[[1]]$transformation <- sprintf("(%s)/(%s)", input[[1]]$transformation, erlang_shape)
			
			# rescale initial value
			if(!is.null(input[[1]]$value)){
				input[[1]]$value <- sprintf("(%s)/(%s)", input[[1]]$value, erlang_shape)
			}

			input_erlang <- rep(input, erlang_shape)
			# remove original input
			inputs[erlang_state] <- NULL

		}

		# set erlang input name
		for(i in seq_along(input_erlang)){
			input_erlang[[i]]$name <- paste0(erlang_state,i)
		}

		inputs <- c(inputs,input_erlang)

	}

	names(inputs) <- NULL

	return(inputs)

}


