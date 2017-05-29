remove_null <- function(x) {

	if(is.list(x)){
		x <- x[!sapply(x,is.null)]
	}

	return(x)
}

clean_args <- function(x) {
	
	x <- x %>% remove_null %>% .[!sapply(.,is.logical) | as.logical(.)] %>% unlist %>% sapply(.,function(xx) str_replace(xx,"TRUE",""))
	
	return(x)
}


get_element <- function(x, key) {

	if(is.list(x)){
		return(sapply(x, function(xx) {xx[[key]]}))
	} else {
		return(x[[key]])
	}

}

get_name <- function(x) {

	get_element(x, "name")

}


find_element <- function(x, key, value) {

	i <- sapply(x, function(xx) {xx[[key]] == value}) %>% which

	return(x[i])

}

#'Export to tracer
#'
#'Export to tracer
#' @param  path character, where to find \code{trace_*.csv}. If \code{NULL} (default), use the \code{path} of the last block (e.g. \code{/pmcmc}).
#' @param  id numeric, indicate which \code{trace_*.csv} to choose. If \code{NULL} (default), use the \code{id} of the last block (default to 0 in SSM).
#' @inheritParams call_ssm
#' @inheritParams plot_X
#' @import readr
#' @export
to_tracer <- function(ssm, path=NULL, id=NULL) {

	if(is.null(path)){

		path <- ssm$hidden$last_path
		
		if(is.null(path)){
			stop("Argument",sQuote("path"),"required", call.=FALSE)	
		}
	}


	if(!is.null(id)){

		df_trace <- sprintf("trace_%s.csv",id) %>% file.path(path,.) %>% read_csv(col_types = cols())

	} else {

		# search for all trace_* in path
		trace_files <- list.files(path) %>% grep("trace_*",.,value=TRUE)

		if(length(trace_files)==0){
			stop("No trace files in directory", dQuote(path))
		}

		if(length(trace_files)>1){
			
			# if more than one, take ssm$summary$id. If missing, send error
			id <- ssm$summary[["id"]]
			if(is.null(id)){
				stop("Use numeric argument",sQuote("id"),"to select one file among:",sQuote(trace_files))
			}
			trace_files <- sprintf("trace_%s.csv",id)

		}

		df_trace <- file.path(path,trace_files) %>% read_csv(col_types = cols())

	}


	df_trace <- data.frame(state=(1:nrow(df_trace))-1,df_trace)
	df_trace$index <- NULL

	# TODO: get id and put it
	write.table(df_trace,file=file.path(path,"tracer.txt"),row.names=FALSE,quote=FALSE,sep="\t")

	invisible(ssm)

}


protect <- function(x) {

	return(paste0("(",x,")"))

}


get_state_variables <- function(reactions) {

	x <- reactions %>% unlist

	x[names(x)%in%c("from","to")] %>% unique %>% return

}


#' Calibrate the number of particles for a SMC
#'
#' This function allows you to optimize the number of particles in your SMC. It displays how the variance and the mean of the likelihood estimator change with increasing number of particles. The idea is then to select an intermediate number of particles such that the estimator is stable, thus saving further computation burden in your PMCMC.
#' @param  n_parts numeric vector containing the number of particles to test.
#' @param  n_replicates for each value of \code{n_parts}, the number of replicates. The higher the more accurate the estimatation of the variance but it will take longer to compute.
#' @param  plot  logical, if \code{TRUE} (default) the function will print 2 diagnostic plots.
#' @param  ... other parameters to be passed to \code{\link{SMC}}.
#' @inheritParams call_ssm
#' @note This function will always run an SMC with PSR approximation since it's the only situation where you will be interested in calibrating the number of particles.
#' @export
#' @import ggplot2 dplyr tidyr
#' @return a list of 3 elements:
#'\itemize{
#'	\item \code{smc} a \code{data_frame} containing one \code{ssm} object per SMC run.
#'	\item \code{summary} a \code{data_frame} containing the SSM summary of each SMC run (loglikelihood etc.).
#'	\item \code{plot} a list of 2 plots: \code{boxplot} and \code{mean_var}
#'}
calibrate_smc <- function(ssm, n_parts, n_replicates, plot = TRUE, ...) {


	df_smc <- crossing(n_parts = unique(n_parts), replicate = 1:n_replicates) %>% group_by(n_parts, replicate) %>% do(ssm = smc(ssm, approx = "psr", n_parts = .$n_parts, ...))

	df_summary <- df_smc %>% group_by(n_parts, replicate) %>% do(as_data_frame(t(.$ssm[[1]]$summary))) %>% ungroup

	df_plot <- df_summary %>% select(-n_data, -n_parameters, -id, -replicate) %>% gather(variable, value, -n_parts)

	p <- ggplot(df_plot, aes(x = n_parts, y = value, group = n_parts)) + facet_wrap(~variable, scales = "free_y")		
	p <- p + geom_boxplot() + geom_jitter()
	p <- p + theme_minimal()
	p_boxplot <- p
	
	if(plot){
		quartz()
		print(p_boxplot)
	}

	df_mean_var <- df_summary %>% select(n_parts, log_likelihood) %>% group_by(n_parts) %>% summarize(variance = var(log_likelihood), mean = mean(log_likelihood))

	df_plot <- df_mean_var %>% gather(variable, value, -n_parts)

	p <- ggplot(df_plot, aes(x = n_parts, y = value)) + facet_wrap(~variable, scales = "free_y")		
	p <- p + geom_line() + geom_point()
	p <- p + theme_minimal()
	p_mean_var <- p

	if(plot){
		quartz()
		print(p_mean_var)
	}


	ans <- list(smc = df_smc, summary = df_summary, plot = list(boxplot = p_boxplot, mean_var = p_mean_var))
	
	invisible(ans)

}



#'Compute hat for X
#'
#'Compute quantile intervals (hat) for states trajectories $X_t$.
#' @param df_X A dataframe with 4 columns: date, index, variable and value. 
#' @param  hat A numeric vector containing the quantiles to compute, centered on the median. For example \code{hat = c(0.5, 0.95)} will compute the interquartile range as well as the 95% credible interval.
#' @export
#' @return A dataframe
get_hat <- function(df_X, hat) {


	prob <- c((1-hat)/2,(1+hat)/2) %>% unique %>% sort 
	dots_summarize <- as.list(sprintf("stats::quantile(value, %s, type=1, na.rm = TRUE)",prob))	
	dots_group_by <- setdiff(names(df_X), c("value","index"))

	hat_label <- paste0(sort(hat)*100,"%")
	dots_names <- c(sprintf("lower_%s",rev(hat_label)),sprintf("upper_%s",hat_label))

	df_hat <- df_X %>% group_by_(.dots=dots_group_by) %>% summarize_(.dots=setNames(dots_summarize,dots_names)) %>% ungroup %>% gather(tmp, value, matches("lower|upper")) %>% separate(tmp,c("hat","level"),sep="_") %>% spread(hat, value)

	return(df_hat)
}


# df_tracer <- as.data.frame(my_mcmc_burn_thin_combined)


# get_args <- function(args_names) {

# 	browser()

# 	arg_names <- formals(fun=sys.function(which=2)) %>% names

# 	env <- parent.frame()

# 	args <- sapply(arg_names, get, envir=env)

# 	return(args)

# }