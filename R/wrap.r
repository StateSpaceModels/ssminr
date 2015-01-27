#'Call to SSM
#'
#'This is the core function to call SSM from R
#' @param ssm a \code{ssm} object, returned by \code{\link{new_ssm}}.
#' @param approx character, approximation used to simulate \code{ssm}:
#' \itemize{
#' 	\item \code{"ode"} deterministic approximation, based on Ordinary Differential Equations.
#' 	\item \code{"sde"} diffusion approximation, based on Stochastic Differential Equations. 
#' 	\item \code{"psr"} Euler-multinomial approximation, based on Poisson process with Stochastic Rates. 
#' }
#' @param do character, algorithm to perform on \code{ssm}
#' @param env calling environment
#' @import dplyr parallel
#' @return a \code{ssm} object updated with latest SSM output and ready to be piped into another SSM block.
call_ssm <- function(ssm, approx=c("ode","sde","psr"), do=c("kalman","kmcmc","ksimplex","mif","pmcmc","simplex","simul","smc"), env=parent.frame()) {

	if(!inherits(ssm,"ssm")){
		stop(sQuote("ssm"),"is not an object of class ssm")
	}

	approx <- match.arg(approx)
	do <- match.arg(do)

	# capture arguments of calling function (remove ssm)
	arg_names <- formals(fun=do) %>% names %>% setdiff(c("ssm","approx"))

	# eval arguments from calling environement
	cmd_args <- sapply(arg_names, function(x) {eval(as.symbol(x), envir=env)})

	if(is.null(cmd_args$root)){
		# default simul directory
		cmd_args$root <- file.path(ssm$model_path, do)
		dir.create(cmd_args$root, showWarnings=FALSE)
	}

	if(!is.null(cmd_args$n_thread) && cmd_args$n_thread=="max"){
		cmd_args$n_thread <- parallel::detectCores()
	}

	theta_in <- r2ssm_resources(ssm$theta, ssm$covmat)

	# write theta_in
	theta_in_path <- file.path(cmd_args$root,sprintf("theta_in_%s.json", cmd_args$id))
	write(rjson::toJSON(theta_in), file=theta_in_path)

	theta_out_path <- file.path(cmd_args$root,sprintf("theta_out_%s.json", cmd_args$id))

	cmd <- sprintf("cd %s/bin; cat %s | ./%s %s %s > %s", ssm$model_path, theta_in_path, do, approx, r2ssm_args(cmd_args), theta_out_path)

	# execute cmd 
	system(cmd)

	# update theta, covmat and summary
	output <- ssm2r_resources(rjson::fromJSON(file=theta_out_path))
	
	ssm$theta <- output$theta
	ssm$covmat <- output$covmat
	if(!is.null(output$summary)){
		ssm$summary <- output$summary		
	}

	invisible(ssm)
}


#'Simulate ssm
#'
#'Function to simulate a \code{ssm}.
#' @param  dt numeric, integration time step.
#' @param  id integer, unique integer identifier that will be appended to the output.
#' @param  root character, root path for output files (if any) (no trailing slash). If \code{NULL} (default), outputs are written in "your_model_path/the_name_of_the_wrapper". 
#' @param  n_thread numeric, number of threads to be used. Default to 1. Use \code{"max"} to set it automatically to the number of cores available on the machine using \code{\link[parallel]{detectCores}}.
#' @param  n_parts numeric, number of particles.
#' @param  start,end character, starting and ending date of model simulation (in \code{YYYY-MM-DD} format). Default to data range.
#' @param  eps_abs_integ numeric, absolute error for adaptive step-size control.
#' @param  eps_rel_integ numeric, relative error for adaptive step-size control.
#' @param  freeze_forcing character, freeze covariates to their value at specified date (in \code{YYYY-MM-DD} format).
#' @param  freq numeric, for simulations outside the data range, print the outputs (and reset incidences to 0 if any) every specified days.
#' @param  interpolator character, gsl interpolator for covariates
#' @param  verbose logical, print logs (verbose). Default to \code{FALSE}.
#' @param  warning logical, print warnings. Default to \code{FALSE}.
#' @param  no_dem_sto logical, turn off demographic stochasticity (if any). Default to \code{FALSE}.
#' @param  no_white_noise logical, turn off white noises (if any). Default to \code{FALSE}.
#' @param  no_diff logical, turn off diffusions (if any). Default to \code{FALSE}.
#' @param  traj logical, print the trajectories. Default to \code{TRUE}.
#' @param  hat logical, print the state estimates. Default to \code{FALSE}.
#' @param  seed_time logical, seed the random number generator with the current time. Default to \code{TRUE}.
#' @inheritParams call_ssm
#' @export
#' @return a \code{ssm} object updated with latest SSM output and ready to be piped into another SSM block.
simul <- function(ssm, approx=c("ode","sde","psr"), dt=NULL, id=0, root=NULL, n_thread=1, n_parts=NULL, start=NULL, end=NULL, eps_abs_integ=NULL, eps_rel_integ=NULL, freeze_forcing=NULL, freq=NULL, interpolator=NULL, verbose=FALSE, warning=FALSE, no_dem_sto=FALSE, no_white_noise=FALSE, no_diff=FALSE, traj=TRUE, hat=FALSE, seed_time=TRUE) {

	approx <- match.arg(approx)

	# run ssm and return updated ssm
	invisible(call_ssm(ssm=ssm, approx=approx, do="simul"))
	
}


#'Run a Kalman filter
#'
#'Function to run a Kalman filter on the \code{sde} approximation of a \code{ssm}.
#' @param  n_obs numeric, number of observations to be fitted (for tempering). If \code{NULL} (default), all observations are fitted.
#' @param  like_min numeric, particles with likelihood smaller than \code{like_min} are considered lost. If \code{NULL} (default) lower bound on likelihood based on machine precision.
#' @param  trace logical, print the trace. Default to \code{TRUE}.
#' @param  diag logical, print the diagnostics outputs (e.g. prediction residuals). Default to \code{TRUE}.
#' @param  prior logical, add log(prior) to the estimated log-likelihood. Default to \code{TRUE}.
#' @inheritParams call_ssm
#' @inheritParams simul
#' @export
#' @return a \code{ssm} object updated with latest SSM output and ready to be piped into another SSM block.
kalman <- function(ssm, dt=NULL, id=0, root=NULL, n_obs=NULL, eps_abs_integ=NULL, eps_rel_integ=NULL, like_min=NULL, freeze_forcing=NULL, interpolator=NULL, verbose=FALSE, warning=FALSE, no_dem_sto=FALSE, no_white_noise=FALSE, no_diff=FALSE, traj=TRUE, hat=FALSE, trace=TRUE, diag=TRUE, prior=FALSE, seed_time=TRUE) {

	# run ssm and return updated ssm
	invisible(call_ssm(ssm=ssm, approx="sde", do="kalman"))
	
}

#'Run a Kalman-MCMC
#'
#'Function to run a MCMC using a Kalman filter to evaluate the log-likelihood. Only for the \code{sde} approximation of a \code{ssm}.
#' @param  iter numeric, number of iterations.
#' @param  cooling numeric, cooling factor (for sampling covariance live tuning or \code{mif} cooling).
#' @param  switch numeric, select switching iteration from initial covariance to empirical one (\code{kmcmc} and \code{pmcmc}) or to update formula introduced in Ionides et al. 2006 (\code{mif})
#' @param  eps_switch numeric, select number of burn-in iterations before tuning epsilon.
#' @param  eps_max numeric, maximum value allowed for epsilon.
#' @param  smooth logical, tune epsilon with the value of the acceptance rate obtained with exponential smoothing. Default to \code{FALSE}.
#' @param  alpha numeric, smoothing factor of exponential smoothing used to compute smoothed acceptance rate (low values increase degree of smoothing)
#' @param  n_traj numeric, number of trajectories stored. If \code{NULL} (default), one trajectory is stored at each iteration.
#' @param  acc logical, print the acceptance rate. Default to \code{TRUE}.
#' @inheritParams call_ssm
#' @inheritParams simul
#' @inheritParams kalman
#' @export
#' @return a \code{ssm} object updated with latest SSM output and ready to be piped into another SSM block.
kmcmc <- function(ssm, dt=NULL, id=0, root=NULL, iter=NULL, n_obs=NULL, cooling=NULL, switch=NULL, eps_switch=NULL, eps_max=NULL, eps_abs_integ=NULL, eps_rel_integ=NULL, smooth=FALSE, alpha=NULL, like_min=NULL, freeze_forcing=NULL, interpolator=NULL, verbose=FALSE, warning=FALSE, no_dem_sto=FALSE, no_white_noise=FALSE, no_diff=FALSE, traj=TRUE, n_traj=NULL, hat=FALSE, trace=TRUE, acc=TRUE, seed_time=TRUE) {

	# run ssm and return updated ssm
	invisible(call_ssm(ssm=ssm, approx="sde", do="kmcmc"))
	
}

#'Run a Kalman-Simplex
#'
#'Function to run a Simplex using a Kalman filter to evaluate the log-likelihood. Only for the \code{sde} approximation of a \code{ssm}.
#' @param size numeric, simplex size used as stopping criteria
#' @export
#' @inheritParams call_ssm
#' @inheritParams simul
#' @inheritParams kalman
#' @inheritParams kmcmc
#' @return a \code{ssm} object updated with latest SSM output and ready to be piped into another SSM block.
ksimplex <- function(ssm, dt=NULL, id=0, root=NULL, iter=NULL, n_obs=NULL, size=NULL, eps_abs_integ=NULL, eps_rel_integ=NULL, like_min=NULL, freeze_forcing=NULL, interpolator=NULL, verbose=FALSE, warning=FALSE, no_dem_sto=FALSE, no_white_noise=FALSE, no_diff=FALSE, trace=TRUE, prior=FALSE, seed_time=TRUE) {

	# run ssm and return updated ssm
	invisible(call_ssm(ssm=ssm, approx="sde", do="ksimplex"))
	
}


#'Run a MIF
#'
#'Function to run a Maximum-likelihood via Iterated Filtering algorithm.
#' @param heat numeric, re-heating across MIF iterations (scales standard deviation of proposals)
#' @param lag numeric, lag for fixed-lag smoothing (proportion of the data)
#' @param ic_only logical, only fit the initial condition using fixed lag smoothing. Default to \code{FALSE}.
#' @export
#' @inheritParams call_ssm
#' @inheritParams simul
#' @inheritParams kalman
#' @inheritParams kmcmc
#' @return a \code{ssm} object updated with latest SSM output and ready to be piped into another SSM block.
mif <- function(ssm, approx=c("ode","sde","psr"), dt=NULL, id=0, root=NULL, n_parts=NULL, iter=NULL, n_obs=NULL, cooling=NULL, switch=NULL, eps_switch=NULL, eps_max=NULL, eps_abs_integ=NULL, eps_rel_integ=NULL, like_min=NULL, heat=NULL, lag=NULL, freeze_forcing=NULL, interpolator=NULL, verbose=FALSE, warning=FALSE, no_dem_sto=FALSE, no_white_noise=FALSE, no_diff=FALSE, traj=TRUE, trace=TRUE, diag=TRUE, prior=FALSE, ic_only=FALSE, seed_time=TRUE) {

	approx <- match.arg(approx)

	# run ssm and return updated ssm
	invisible(call_ssm(ssm=ssm, approx=approx, do="mif"))
	
}



#'Run a Particle-MCMC
#'
#'Function to run a MCMC using a SMC (particle filter) to evaluate the log-likelihood.
#' @inheritParams call_ssm
#' @inheritParams simul
#' @inheritParams kalman
#' @inheritParams kmcmc
#' @export
#' @return a \code{ssm} object updated with latest SSM output and ready to be piped into another SSM block.
pmcmc <- function(ssm, approx=c("ode","sde","psr"), dt=NULL, id=0, root=NULL, iter=NULL, n_parts=NULL, n_thread=1, n_obs=NULL, cooling=NULL, switch=NULL, eps_switch=NULL, eps_max=NULL, eps_abs_integ=NULL, eps_rel_integ=NULL, smooth=FALSE, alpha=NULL, like_min=NULL, freeze_forcing=NULL, interpolator=NULL, verbose=FALSE, warning=FALSE, no_dem_sto=FALSE, no_white_noise=FALSE, no_diff=FALSE, traj=TRUE, n_traj=NULL, hat=FALSE, trace=TRUE, acc=TRUE, seed_time=TRUE) {

	approx <- match.arg(approx)

	# run ssm and return updated ssm
	invisible(call_ssm(ssm=ssm, approx=approx, do="pmcmc"))
	
}


#'Run a Simplex
#'
#'Function to run a Simplex on the \code{ode} approximation of a \code{ssm}.
#' @param least_squares logical, minimize the sum of squared errors instead of maximizing the likelihood. Default to \code{FALSE}.
#' @export
#' @inheritParams call_ssm
#' @inheritParams simul
#' @inheritParams kalman
#' @inheritParams ksimplex
#' @inheritParams kmcmc
#' @return a \code{ssm} object updated with latest SSM output and ready to be piped into another SSM block.
simplex <- function(ssm, dt=NULL, id=0, root=NULL, iter=NULL, n_obs=NULL, size=NULL, eps_abs_integ=NULL, eps_rel_integ=NULL, like_min=NULL, freeze_forcing=NULL, interpolator=NULL, verbose=FALSE, warning=FALSE, no_dem_sto=FALSE, no_white_noise=FALSE, no_diff=FALSE, trace=TRUE, prior=FALSE, least_squares=FALSE, seed_time=TRUE) {

	# run ssm and return updated ssm
	invisible(call_ssm(ssm=ssm, approx="ode", do="simplex"))
	
}



#'Run a SMC
#'
#'Function to run a Sequential Monte-Carlo algorithm on a \code{ssm}
#' @param no_filter logical, do not filter. Default to \code{FALSE}.
#' @export
#' @inheritParams call_ssm
#' @inheritParams simul
#' @inheritParams kalman
#' @inheritParams kmcmc
#' @return a \code{ssm} object updated with latest SSM output and ready to be piped into another SSM block.
smc <- function(ssm, approx=c("ode","sde","psr"), dt=NULL, id=0, root=NULL, n_parts=NULL, n_thread=1, iter=NULL, n_obs=NULL, eps_max=NULL, eps_abs_integ=NULL, eps_rel_integ=NULL, like_min=NULL, freeze_forcing=NULL, interpolator=NULL, verbose=FALSE, warning=FALSE, no_dem_sto=FALSE, no_white_noise=FALSE, no_diff=FALSE, traj=TRUE, hat=FALSE, trace=TRUE, diag=TRUE, prior=FALSE, no_filter=FALSE, seed_time=TRUE) {

	approx <- match.arg(approx)

	# run ssm and return updated ssm
	invisible(call_ssm(ssm=ssm, approx=approx, do="smc"))
	
}


# generate_lhs <- function(sample_size=100) {

# 	dir_lhs <- file.path(dir_model,"lhs")

# 	if(!file.exists(dir_lhs)){
# 		dir.create(dir_lhs)
# 	}

# 	# remove dirac prior
# 	ssm_model <- readRDS(file.path(dir_model,"ssm_model.rds"))


# 	theta_estimated_names <- names(ssm_model$theta$init)
# 	covmat <- ssm_model$theta$covmat
# 	theta_priors <- ssm_model$theta_priors

# 	for(i in 1:sample_size){

# 		init_theta <- sample_from_prior(theta_priors, theta_estimated_names)
# 		names(init_theta) <- theta_estimated_names
# 		# init_theta[["beta_I"]] <- runif(1,0,2)
# 		# init_theta[["beta_H"]] <- runif(1,1,5)
# 		# init_theta[["beta_F"]] <- runif(1,0,2)
# 		# init_theta[["shift_sigmo_H"]] <- runif(1,35,45)
# 		# init_theta[["shape_sigmo_H"]] <- runif(1,0,1)
# 		# init_theta[["shift_sigmo_com"]] <- runif(1,45,55)
# 		# init_theta[["shape_sigmo_com"]] <- runif(1,0,1)

# 		theta_json <- list(resources=ssm_resources(init_theta, covmat))

# 		write(toJSON(theta_json),file=file.path(dir_lhs,paste0("theta_",i-1,".json")))

# 	}
# }

# do_lhs <- function(n_cores=NULL) {

# 	dir_lhs <- file.path(dir_model,"lhs")
# 	dir_bin <- file.path(dir_model,"bin")

# 	# look into dir_model/lhs and extract all theta.json
# 	theta_files <- grep("theta_[0-9]+.json",list.files(dir_lhs),value=TRUE)

# 	# # run a for loop with system call
# 	for(i in seq_along(theta_files)){

# 		theta_file <- theta_files[i]
# 		cmd <- sprintf("cd %s; cat %s | ./simplex -M 10000 --prior > %s/theta_map_simplex_%s.json", 
# 			dir_bin, 
# 			file.path(dir_lhs,theta_file),
# 			dir_lhs,
# 			i-1
# 			)

# 		system(cmd, wait=FALSE)

# 	}

# }


# summarize_lhs <- function() {

# 	dir_lhs <- file.path(dir_model,"lhs")
# 	dir_rds <- file.path(dir_lhs,"rds")

# 	if(!file.exists(dir_rds)){
# 		dir.create(dir_rds)
# 	}

# 	file_names <- grep(".*simplex.*",list.files(dir_lhs),value=TRUE)

# 	# read all files
# 	names(file_names) <- file_names
# 	df_ll <- ldply(file_names, function(file_name) {
# 		res <- fromJSON(file=file.path(dir_lhs,file_name))
# 		return(res$resources[[3]]$data$log_ltp)
# 	}, .progress="text")

# 	names(df_ll) <- c("file_name","log_ltp")

# 	df_theta <- ddply(df_ll,"file_name",function(df) {

# 		res <- fromJSON(file=file.path(dir_lhs,df$file_name))	
# 		unlist(res$resources[[1]]$data)

# 	}, .progress="text")

# 	df_summarize <- left_join(df_theta, df_ll, by="file_name")

# 	saveRDS(df_summarize, file=file.path(dir_rds, "summarize_lhs.rds"))

# 	if(0){
# 		# analysis
# 		x <- readRDS(file.path(dir_model,"lhs_simplex","rds","summarize_lhs.rds"))
# 		x <- df_summarize
# 		min_log_ltp <- -2180
# 		df_plot <- x %>% filter(log_ltp > min_log_ltp & log_ltp < -130) %>% gather(variable, value, -file_name)
# 		p <- ggplot(df_plot, aes(x=value)) + facet_wrap(~variable, scales="free")
# 		p <- p + geom_histogram()
# 		print(p)

# 		# df_plot <- x %>% filter(log_ltp > min_log_ltp)
# 		# p <- ggplot(df_plot, aes(x=R0__sierra_leone, y=t_0__sierra_leone))+geom_point()
# 		# print(p)

# 		dir <- file.path(dir_model,"lhs_simplex","figures")
# 		if(!file.exists(dir)){
# 			dir.create(dir, recursive=TRUE)
# 		}
# 		ggsave(file.path(dir,paste0("dist_simplex_ll>",min_log_ltp,".pdf")))
# 	}

# }

# get_map_lhs <- function() {

# 	dir_lhs <- file.path(dir_model,"lhs")
# 	dir_rds <- file.path(dir_lhs,"rds")

# 	ssm_model <- readRDS(file.path(dir_model,"ssm_model.rds"))

# 	df_lhs <- readRDS(file.path(dir_rds,"summarize_lhs.rds")) 

# 	if(all(df_lhs$log_ltp <= 0)){
# 		df_lhs <- df_lhs %>% filter(log_ltp < -10)
# 	}

# 	df_lhs <- df_lhs %>% filter(log_ltp==max(log_ltp))

# 	cmd <- sprintf("cp %s %s/theta_map_simplex.json",file.path(dir_lhs, df_lhs$file),dir_model)
# 	system(cmd)

# }

# do_mcmc <- function(n_iter_short_run=1e4, n_iter_long_run=1e5, n_traj=1000) {

# 	dir_bin <- file.path(dir_model,"bin")
# 	dir_mcmc <- file.path(dir_model,"mcmc")
# 	if(!file.exists(dir_mcmc)){
# 		dir.create(dir_mcmc)
# 	}

# 	cmd <- sprintf("cd %s; cat %s/theta_map_simplex.json | ./pmcmc --iter %0.f --eps_switch 20 --cooling 0.99 --switch 500 --seed_time | ./pmcmc --iter %0.f --eps_switch 20 --cooling 0.99 --switch 500 --seed_time --trace --traj --n_traj %s --acc --root %s > %s/theta_mean_mcmc.json
# 		",
# 		dir_bin,
# 		dir_model,
# 		n_iter_short_run,
# 		n_iter_long_run,
# 		n_traj,
# 		dir_mcmc,
# 		dir_mcmc
# 		)

# 	system(cmd)

# }

