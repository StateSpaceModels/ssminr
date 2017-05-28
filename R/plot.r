
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
plot_X <- function(ssm, path=NULL, id=NULL, stat=c("none", "median", "mean"), hat=NULL, scales="free_y", fit_only=FALSE, collapse_erlang=TRUE) {

	stat <- match.arg(stat)

	data_colour <- "#1f78b4"

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

	obs_var <- get_name(ssm$observations)
	ran_obs_var <- sprintf("ran_%s", obs_var)

	# TODO find accumulators and bind them to state_variables and search in obs parameter (mean etc)
	# accumulators should be extracted at compilation 
	# here for simplicity we assume that observation names are of the form state_obs so state can easily be retrieved
	obs_state <- obs_var[str_detect(obs_var, "_obs")] %>% str_replace("_obs", "")

	obs_dist <- get_element(ssm$observations, "distribution")
	names(obs_dist) <- obs_var

	# keep observed states and observations at the data dates
	
	df_data_obs <- ssm$data %>% rename(state = time_series)
	df_data_ran_obs <- df_data_obs %>% mutate(state = sprintf("ran_%s", state))
	df_data_state <- df_data_obs %>% filter(str_detect(state, "_obs")) %>% mutate(state = str_replace(state, "_obs", ""))

	df_data_date <- bind_rows(df_data_obs, df_data_ran_obs, df_data_state)

	df_X_obs <- df_X %>% filter(state %in% c(obs_var, obs_state, ran_obs_var)) %>% semi_join(df_data_date, c("date", "state"))

	df_X <- df_X %>% anti_join(df_X_obs, "state") %>% bind_rows(df_X_obs)


	
# remove binomial observations
	if(any(obs_dist == "binomial")) {

		obs_binomial <- obs_dist[obs_dist == "binomial"]

		df_X_binomial <- df_X %>% filter(state %in% c(names(obs_binomial), sprintf("ran_%s", names(obs_binomial))))

		df_X <- df_X %>% anti_join(df_X_binomial, c("state"))

		# binomial denominator
		n_tested_name <- get_element(ssm$observation, "n")
		names(n_tested_name) <- obs_var
		n_tested_name <- unlist(n_tested_name)


		df_n <- data_frame(n_name = n_tested_name, state = names(n_tested_name)) %>% 
		group_by(n_name, state) %>% 
		do(n = find_element(ssm$inputs, "name", .$n_name)[[1]]$value, date = find_element(ssm$inputs, "name", .$n_name)[[1]]$value %>% names %>% as.Date) %>% 
		ungroup %>% 
		unnest(n, date) %>% 
		filter(n != 0)

		df_n_ran <- df_n %>% mutate(state = sprintf("ran_%s", state))

		df_n <- df_n %>% bind_rows(df_n_ran)

		
		# binomial posterior
		df_prop <- df_X_binomial %>% left_join(df_n, by = c("date", "state")) %>% mutate(p = value/n*100)

		# binomial data + 95% CI
		df_data_prop <- ssm$data %>% rename(state = time_series) %>% inner_join(df_n, by = c("date", "state")) %>% 
		group_by(date, state) %>% 
		do(binom.confint(.$value, .$n, methods = "exact")) %>% 
		ungroup %>% 
		gather(variable, value, c(mean, lower, upper)) %>% 
		mutate(value  =value*100) %>% 
		spread(variable, value)
		
		## add ran_
		df_data_prop_ran <- df_data_prop %>% mutate(state = sprintf("ran_%s", state))
		df_data_prop <- df_data_prop %>% bind_rows(df_data_prop_ran)

		# plot
		p <- ggplot(df_prop) + facet_wrap(~state, scales = "free_x")
		p <- p + geom_violin(data = df_prop, aes(x = date, y = p, group = date))
		p <- p + geom_pointrange(data =  df_data_prop, aes(x = date, y = mean, ymin = lower, ymax = upper), col = data_colour)
		p <- p + theme_minimal() + ylab("Proprotion (%)") + xlab("Time")

		print(p)
		quartz()

		ssm$plot$X_binom <- p

	}


	# remove states with NA values
	if(any(is.na(df_X$value))){
		
		warnings("NA values in states are removed")

		df_remove_index <- df_X %>% filter(is.na(value))
		df_X <- df_X %>% anti_join(df_remove_index, "index")

	}

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

		df_hat <- get_hat(df_X, hat)
		
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
	
	# plot data on expected and ran observations
	df_data_ran <- df_data %>% mutate(state = sprintf("ran_%s", state))
	df_data <- df_data %>% bind_rows(df_data_ran)

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

		hat_label <- sort(unique(df_hat$level))
		alpha_values <- seq(0.2,0.6,len=length(hat_label)) %>% rev
		names(alpha_values) <- hat_label
		p <- p + geom_ribbon(aes(ymin=lower, ymax=upper, alpha=level)) + scale_alpha_manual("Level", values=alpha_values)
	}

	if(stat!="none"){
		# add stat
		p <- p + geom_line(data=df_stat, aes(y=value))
	}

	p <- p + geom_point(data=df_data, aes(y=value), shape = 1, colour = data_colour)

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


#'Plot priors
#'
#'Plot prior distribution as specified in \code{ssm} object.
#' @param theta_names character, specify which theta priors to plot. By default (=\colde{NULL}) all theta priors are plotted.
#' @param quantile_limits numeric, vector of length 2 specifying the limits of the priors in terms of quantile. Default to the 1% and 99% quantiles (=\code{c(0.01, 0.99)}).
#' @param x_limits numeric, vector of length 2 specifying the limits of the priors in terms of theta. Default to \code{NULL}.
#' @param plot logical, if \code{FALSE} a dataframe of prior values will be returned insread iof a plot
#' @inheritParams call_ssm
#' @inheritParams plot_X
#' @export
#' @import dplyr ggplot2
#' @return either a \code{ggplot} or a \code{data.frame}
plot_priors <- function(ssm, theta_names=NULL, quantile_limits=c(0.01,0.99), x_limits=NULL, plot=TRUE){

	priors <- ssm$priors
	prior_names <- sapply(priors, function(x) {x$name})
	names(priors) <- prior_names
	
	if(!is.null(theta_names)){
		select_priors <- intersect(theta_names,names(priors))
		priors <- priors[select_priors]
	}

	df_prior_names <- data_frame(theta=prior_names)

	compute_plot <- function(prior_name) {

		prior <- priors[[prior_name]]
		x_range <- x_limits[[prior_name]]

		if(is.null(x_range)){
			x_range <- rep(NA,2)
		}

		if(any(is.na(x_range))){
		# use quantile limits
			qrior_dist <- paste0("q",prior$dist)
			x_range[is.na(x_range)] <- do.call(qrior_dist,args=c(prior$args,list(p=quantile_limits[is.na(x_range)])))
		}	

		x <- seq(min(x_range),max(x_range),length=1000)

		dprior <- paste0("d",prior$dist)
		density <- do.call(dprior,c(prior$args,list(x=x)))
		
		return(data_frame(x=x, density=density))
	}

	df_dist <- group_by(df_prior_names, theta) %>% do(compute_plot(unlist(.))) %>% ungroup

	if(plot){
		p <- ggplot(df_dist,aes(x=x,y=density))+facet_wrap(~theta,scales="free")+geom_line()
		print(p + theme_minimal())		
	} else {
		return(df_dist)
	}
}


#'Plot priors and posteriors
#'
#'Plot priors and posteriors distribution for all estimated parameters.
#' @inheritParams call_ssm
#' @inheritParams plot_X
#' @export
#' @import ggplot2 tidyr dplyr readr
#' @return a \code{ggplot} object
plot_theta <- function(ssm, path=NULL, id=NULL, x_limits = NULL) {

	if(is.null(path)){

		path <- ssm$hidden$last_path

		if(is.null(path)){
			stop("Argument",sQuote("path"),"required", call.=FALSE)	
		}
	}

	if(!is.null(id)){

		trace_file <- sprintf("trace_%s.csv",id)

	} else {

		# search for all X_* in path
		trace_file <- list.files(path) %>% grep("trace_[0-9]+.csv",.,value=TRUE)

		if(length(trace_file)==0){
			stop("No trace files in directory", dQuote(path))
		}

		if(length(trace_file)>1){

			# if more than one, take ssm$summary$id. If missing, send error
			id <- ssm$summary[["id"]]
			if(is.null(id)){
				stop("Use numeric argument",sQuote("id"),"to select one file among:",sQuote(trace_file))
			}
			trace_file <- sprintf("trace_%s.csv",id)

		}

	}


	df_posterior <- read_csv(file.path(path, trace_file), col_types = cols(.default = col_guess())) %>% 
	gather(theta, value, -index, -fitness) %>% mutate(distribution = "posterior")  

	# df_limits <- df_posterior %>% group_by(theta) %>% summarize(theta_min = min(value), theta_max = max(value))
	
	# compute prior
	df_prior <- plot_priors(ssm, quantile_limits=c(0.0001,0.9999), x_limits = x_limits, plot=FALSE) %>% mutate(distribution = "prior")  
	
	p <- ggplot(df_prior, aes(fill = distribution)) + facet_wrap(~theta, scales="free")
	p <- p + geom_area(data=df_prior, aes(x=x, y=density), alpha=0.4)
	p <- p + geom_histogram(data=df_posterior, aes(x=value, y=..density..),alpha=0.6)
	p <- p + scale_fill_discrete("Distribution", breaks=c("prior","posterior"))
	p <- p + xlab("Theta") + ylab("Density")
	p <- p + theme_minimal() + theme(legend.position="top")
	print(p)

	# add to ssm plot
	ssm$plot$theta <- p

	invisible(ssm)
	
}


