
#'Interactive display of SSM model
#'
#'This function displays the SSM model as a force network or a diagramme in a web browser, allowing for manipulation.
#' @param ssm a \code{\link{ssm}} object
#' @param collapse_erlang logical, if \code{TRUE} erlang compartments are collapsed into a single compartment to improve visibility.
#' @param display logical, if \code{TRUE} erlang compartments are collapsed into a single compartment to improve visibility.
#' @param direction character, direction of the flowchart (\code{display=="diagramme"} only). Available options are \code{LR} (left to right, by default), \code{TB} (top to bottom), \code{RL} and \code{BT}.
#' @export
#' @import networkD3 DiagrammeR
plot_model <- function(ssm, collapse_erlang = TRUE, display=c("diagramme", "network"), engine = "dot", label = TRUE, label_size_max=Inf, simplify_state_names = TRUE, direction = "LR") {

	display <- match.arg(display)

	from <- get_element(ssm$reactions, "from")
	to <- get_element(ssm$reactions, "to")
	rate <- get_element(ssm$reactions, "rate")
	network_data <- data_frame(reaction= seq_along(from), from, to, rate)

	if(collapse_erlang & !is.null(ssm$erlang_shapes)){

		erlang_shapes <- ssm$erlang_shapes
		erlang_states <- names(erlang_shapes)
		df_erlang <- data_frame(from= erlang_states, shape=erlang_shapes)

		remove_to <- sapply(erlang_states, function(erlang_state) erlang_state %>% erlang_name(2:(erlang_shapes[erlang_state]))) %>% unlist

		revalue_from <- erlang_name(erlang_states, 1)
		names(revalue_from) <- sapply(erlang_states, function(erlang_state) erlang_state %>% erlang_name((erlang_shapes[erlang_state]))) %>% unlist

		revalue_sum <- erlang_states
		names(revalue_sum) <- sapply(erlang_states, function(erlang_state) erlang_state %>% erlang_name(1:(erlang_shapes[erlang_state])) %>% paste(collapse=" + ")) %>% unlist

		network_data <- network_data %>% mutate(from=revalue(from, revalue_from)) %>%
		filter(!to%in%remove_to) %>% gather(type, state, -c(reaction, rate)) %>% 
		mutate(state=str_replace_all(state, "__erlang_[0-9]+", "")) %>% spread(type, state) 

		# replace sum of erlang compartments in the rates
		rate <- network_data$rate
		for(sum_name in names(revalue_sum)){
			# need to use fixed() because of the special charcaters
			rate <- str_replace_all(rate, fixed(sum_name), revalue_sum[[sum_name]]) 
		}
		network_data$rate <- rate	
		
		# if from is an erlang compartment => divide rate 
		network_data <- network_data %>% left_join(df_erlang, by="from") %>% replace_na(list(shape=1)) %>% mutate(rate = sprintf("(%s)/%s", rate, shape)) %>% select(-shape)

		# simplify	
		all_var <- c(get_name(my_ssm$inputs),revalue_sum %>% unname)
		network_data <- network_data %>% mutate(rate = sympy_simplify(rate, all_var))		

	}


	if(display=="diagramme"){

		df_nodes <- network_data %>% .[c("from","to")] %>% unlist %>% unique %>% data_frame(node=., state=node, original=node)

		if(simplify_state_names){
			# improve visibility

			# simplify state names

			if(any(str_detect(df_nodes$state, erlang_name()))){

				#simplify erlang
				df_erlang <- df_nodes %>% filter(str_detect(state, erlang_name())) %>% separate(state, c("state", "erlang"), sep=erlang_name())
				df_nodes <- df_nodes %>% filter(!str_detect(state, erlang_name())) %>% bind_rows(df_erlang)
			} else {
				df_nodes$erlang <- NA
			}

			if(any(str_detect(df_nodes$state, pop_name()))){

				#simplify pop
				df_pop <- df_nodes %>% filter(str_detect(state, pop_name())) %>% separate(state, c("state", "pop"), sep=pop_name())
				df_nodes <- df_nodes %>% filter(!str_detect(state, pop_name())) %>% bind_rows(df_pop)

			} else {

				df_nodes$pop <- NA
			}			

			df_nodes <- df_nodes %>% group_by(state, pop, erlang) %>% 
			mutate(
				sub = paste(na.omit(c(pop, erlang)), collapse=","),
				label = sprintf("%s<SUB>%s</SUB>", state, sub) %>% str_replace_all("<SUB></SUB>",""), 
				node = sprintf("%s [label = <%s>]", original, label)
				) 

			## replace in rates
			replace_state <- df_nodes$label
			names(replace_state) <- sprintf("\\b%s\\b",df_nodes$original)

			network_data <- network_data %>% mutate(rate = rate %>% str_replace_all(replace_state))

		} 

		node_statement <- df_nodes$node %>% paste(collapse="; ")

		if(is.finite(label_size_max)){
			network_data <- network_data %>% mutate(rate = rate %>% str_sub(start=1L, end=label_size_max) %>% paste0("..."))
		}

		if(!label){

			network_data <- network_data %>% mutate( rate = "")

		}

		edge_statement <- network_data %>% unite(edge, c(from, to), sep="->") %>% 
		mutate(edge = sprintf("%s [label = <%s>]", edge, rate)) %>% 
		select(edge) %>% unlist %>% unname %>% paste(collapse=" ") %>% convert_symbol(to="html")			
		
		gviz_cmd <- sprintf("
			digraph circles {

 		 # a 'graph' statement
				graph [overlap = true, fontsize = 10, rankdir = %s]

 		 # several 'node' statements
				node [shape = circle,
				fontname = Helvetica,
				style = filled,
				color = grey,
				fillcolor = lightsteelblue]
				%s

 		 # several 'edge' statements
				edge[color = grey]
				%s
			}
			", direction, node_statement, edge_statement)

		return(grViz(gviz_cmd, engine = engine))

	}


	if(display=="network"){

		return(simpleNetwork(network_data, Source="from", Target="to", zoom = TRUE))
	}

}


#'Plot states
#'
#'Function to plot state trajectories contained in the \code{X_*.csv} files generated by SSM. Pipeable.
#' @param  path character, where to find \code{X_*.csv}. If \code{NULL} (default), use the \code{path} of the last block (e.g. \code{/pmcmc}).
#' @param  id numeric, indicate which \code{X_*.csv} to choose. If \code{NULL} (default), use the \code{id} of the last block (default to 0 in SSM).
#' @param  stat character, whether to plot a summary statistics of the state. Either \code{"mean"} or \code{"median"}. Default to \code{"none"}.
#' @param  hat numeric, vector of credible intervals, between 0 and 1, e.g. \code{hat=c(0.5, 0.95)} for 50 and 95\% credible intervals.
#' @param  scales character, should scales be \code{"fixed"}, \code{"free"}, or free in one dimension: \code{"free_x"}, \code{"free_y"} (the default).
#' @param  fit_only logical, whether to show only the fit to the data.
#' @inheritParams call_ssm
#' @export
#' @import ggplot2 tidyr dplyr readr
#' @seealso \code{\link{plot_theta}}
#' @return a \code{ssm} object updated with latest SSM output and ready to be piped into another SSM block.
plot_X <- function(ssm, path=NULL, id=NULL, stat=c("median", "mean", "none"), hat=NULL, scales="free_y", fit_only=FALSE, collapse_erlang=TRUE) {

	stat <- match.arg(stat)

	if(is.null(path)){

		path <- ssm$hidden$last_path

		if(is.null(path)){
			stop("Argument",sQuote("path"),"required", call.=FALSE)	
		}
	}

	if(!is.null(id)){

		X_file <- sprintf("X_%s.csv",id)

	} else {

		# search for all X_* in path
		X_file <- list.files(path) %>% grep("X_*",.,value=TRUE)

		if(length(X_file)==0){
			stop("No X files in directory", dQuote(path),"..... The Truth is Out There")
		}

		if(length(X_file)>1){

			# if more than one, take ssm$summary$id. If missing, send error
			id <- ssm$summary[["id"]]
			if(is.null(id)){
				stop("Use numeric argument",sQuote("id"),"to select one file among:",sQuote(X_file))
			}
			X_file <- sprintf("X_%s.csv",id)

		}

	}


	# make everything double except date
	# get col names
	col_names <- file.path(path,X_file) %>% read_csv(n_max = 0, col_types = cols()) %>% names
	col_types <- rep("d", length(col_names))
	names(col_types) <- col_names
	col_types["date"] <- "D"
	col_types <- col_types %>% paste(collapse="")

	df_X <- read_csv(file.path(path,X_file), col_types = col_types) %>%  
	mutate(date=as.Date(date)) %>% 
	gather(state, value, -date, -index)  

	if(collapse_erlang){

		df_X <- collapse_erlang(df_X)

	}

	# separate pop only if collapse_erlang. Otherwise we loose erlang order (always last).
	if(collapse_erlang && any(str_detect(df_X$state, pop_name()))){

		df_X_pop <- df_X %>% filter(str_detect(state, pop_name())) %>% separate(state, c("state","pop"), sep=pop_name())
		df_X <- df_X %>% filter(!str_detect(state, pop_name())) %>% bind_rows(df_X_pop) %>% arrange(index, date, pop, state)

	}

	# any stat?
	if(stat!="none"){

		stat <- ifelse(stat=="median","stats::median",stat)
		dots_summarize <- list(sprintf("%s(value)",stat))
		dots_group_by <- setdiff(names(df_X), c("value","index"))
		df_stat <- df_X %>% group_by_(.dots=dots_group_by) %>% summarize_(.dots=setNames(dots_summarize,"value"))
	}

	# any hat?
	if(!is.null(hat)){

		prob <- c((1-hat)/2,(1+hat)/2) %>% unique %>% sort 
		dots_summarize <- as.list(sprintf("stats::quantile(value, %s, type=1, na.rm = TRUE)",prob))	
		dots_group_by <- setdiff(names(df_X), c("value","index"))

		hat_label <- paste0(sort(hat)*100,"%")
		dots_names <- c(sprintf("lower_%s",rev(hat_label)),sprintf("upper_%s",hat_label))

		df_hat <- df_X %>% group_by_(.dots=dots_group_by) %>% summarize_(.dots=setNames(dots_summarize,dots_names)) %>% ungroup %>% gather(tmp, value, matches("lower|upper")) %>% separate(tmp,c("hat","level"),sep="_") %>% spread(hat, value)

	}

	if(!is.null(hat)){
		df_plot <- df_hat
	} else {
		df_plot <- df_X
	}

	df_data <- ssm$data %>% dplyr::rename(state=time_series)

	if(collapse_erlang && any(str_detect(df_data$state, pop_name()))){

		df_data_pop <- df_data %>% filter(str_detect(state, pop_name())) %>% separate(state, c("state","pop"), sep=pop_name())
		df_data <- df_data %>% filter(!str_detect(state, pop_name())) %>% bind_rows(df_data_pop) %>% arrange(date, pop, state)

	}

	by_names <- intersect(names(df_data), names(df_plot)) %>% setdiff("value")
	df_data <- df_data %>% semi_join(df_plot, by=by_names)

	if(fit_only){

		# keep states that match data
		df_plot <- df_plot %>% semi_join(df_data, by_names)
		if(stat!="none"){
			df_stat <- df_stat %>% semi_join(df_data, by_names)			
		}
	}

	facet_formula <- ifelse("pop"%in%names(df_plot), "pop~state", "~state") %>% as.formula
	p <- ggplot(data=df_plot, aes(x=date)) + facet_wrap(facet_formula, scales=scales)



	if(is.null(hat)){
		# plot traj

		# choose alpha
		n_index <- n_distinct(df_plot$index)
		alpha <- ifelse(n_index > 10, min(c(0.1,10/n_index)), 1)
		p <- p + geom_line(aes(y=value, group=index), alpha=alpha)		
	} else {

		alpha_values <- seq(0.2,0.6,len=length(hat_label)) %>% rev
		names(alpha_values) <- hat_label
		p <- p + geom_ribbon(aes(ymin=lower, ymax=upper, alpha=level)) + scale_alpha_manual("Level", values=alpha_values)
	}

	if(stat!="none"){
		# add stat
		p <- p + geom_line(data=df_stat, aes(y=value))
	}

	p <- p + geom_point(data=df_data, aes(y=value))

	print(p + theme_minimal())

	# add to ssm plot
	ssm$plot$X <- p

	invisible(ssm)
}


#'Plot data
#'
#'Plot the data of your \code{ssm} object.
#' @inheritParams call_ssm
#' @inheritParams plot_X
#' @export
#' @import ggplot2 tidyr
plot_data <- function(ssm, scales="free_y") {

	if(!inherits(ssm,"ssm")){
		stop(sQuote("ssm"),"is not an object of class ssm")
	}

	df_data <- ssm$data %>% dplyr::rename(state=time_series)

	if(any(str_detect(df_data$state, pop_name()))){

		df_data <- df_data %>% separate(state, c("state","pop"), sep=pop_name())
		# Not sure in what situation I would have some time-series with _pop_ and some other without..
		# df_data_pop <- df_data %>% filter(str_detect(state, pop_name())) %>% separate(state, c("state","pop"), sep=pop_name())
		# df_data <- df_data %>% filter(!str_detect(state, pop_name())) %>% mutate(pop = ) %>% bind_rows(df_data_pop) %>% arrange(date, pop, state)

	} else {

		df_data <- df_data %>% mutate(pop = ssm$pop)

	}

	p <- ggplot(df_data, aes(x=date, y=value)) + facet_wrap(pop~state, scales=scales)
	p <- p + geom_line() + geom_point() + theme_minimal()
	print(p)

	# add to ssm plot
	ssm$plot$data <- p

	invisible(ssm)

}

# plot_theta <- function(ssm) {

# 	# posterior vs prior distribution of parameters

# 	# get the root of the preceding function



# 	# pass ssm to the next
# 	invisible(ssm)

# }


