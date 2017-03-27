get_ssm_log_like <- function(ssm) {
	
	x <- ssm$summary

	data_frame(log_like=x[str_detect(names(x),"log_*")], ssm=list(ssm))

}

#'Get max LHS
#'
#'Extract the \code{ssm} object with highest log-likelihood from \code{lhs}.
#' @param lhs a \code{tbl}, as returned by \code{\link{do_lhs}}
#' @export
#' @import dplyr
#' @seealso \code{\link{do_lhs}}
#' @return \code{ssm} object
get_max_lhs <- function(lhs) {

	lhs %>% do(get_ssm_log_like(.$ssm)) %>% ungroup %>% filter(log_like==max(log_like)) %>% select(ssm) %>% .[[1,1]]

}

#'Perform LHS
#'
#'Function to perform any SSM block on a LHS. The LHS is currently a random sample from the prior distribution.
#' @param  ... additional arguments to be passed to \code{do}
#' @inheritParams call_ssm
#' @inheritParams sample_prior
#' @export
#' @import dplyr
#' @seealso \code{\link{get_max_lhs}}
#' @return a \code{tbl} containing all \code{ssm} outputed by \code{do}.
do_lhs <- function(ssm, n, do=c("kalman","kmcmc","ksimplex","mif","pmcmc","simplex","simul","smc"), ...) {

	do <- match.arg(do)

	# create lhs directory
	dir_lhs <- ssm$model_path %>% file.path("lhs")

	dir.create(dir_lhs, showWarnings=FALSE)

	# run LHS
	ssm %>% sample_prior(n, method="random") %>% dplyr::do(ssm=do.call(do, list(ssm=.$ssm, id=.$id, root=dir_lhs, ...)))

}


