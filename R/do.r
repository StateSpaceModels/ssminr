


simul <- function(pipe=NULL, model=NULL, theta=NULL, dt=NULL, id=0, root=NULL, suffix=NULL, n_thread=NULL, n_parts=NULL, start=NULL, end=NULL, eps_abs_integ=NULL, eps_rel_integ=NULL, freeze_forcing=NULL, freq=NULL, interpolator=NULL, verbose=FALSE, warning=FALSE, no_dem_sto=FALSE, no_white_noise=FALSE, no_diff=FALSE, traj=TRUE, hat=FALSE, seed_time=TRUE) {

	if(!is.null(pipe)){
		pipe <- eval(pipe)
		theta <- pipe$theta
		model <- pipe$model
	}

	if(is.null(root)){
		# default simul directory
		root <- file.path(model, "simul")
		dir.create(root, showWarnings=FALSE)
	}
	
	if(is.null(theta)){
		# default resource
		theta <- file.path(model, "theta.json")
	}

	# generate cmd
	arg_names <- formals() %>% names %>% setdiff(c("model","theta","pipe"))
	cmd_args <- sapply(arg_names, function(x) {eval(as.symbol(x))})

	cmd <- sprintf("cd %s/bin; cat %s | ./simul %s", model, theta, r2ssm_args(cmd_args))

	# execute cmd
	x <- system(cmd, intern=TRUE)
	
	theta_out <- file.path(root,sprintf("theta_simul_%s.json", id))

	write(x, file=theta_out)

	# return pipe
	invisible(list(theta=theta_out, model=model))
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

