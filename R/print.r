add_line <- function() {

	cat("----------------------\n")

}

print_prior <- function(prior) {

	cat(prior$name,"~",paste0(prior$dist,"("),paste(prior$args %>% names, prior$args, sep=" = ", collapse=", "),")","\n")

}


#'Print ssm
#'
#'Function to print basic information on the \code{ssm}
#' @param x a \code{ssm} object, returned by \code{\link{new_ssm}}.
#' @export
#' @import dplyr
#' @importFrom plyr l_ply
print.ssm <- function(x, ...) {

	ssm <- x

	# model path
	cat("path: ",dQuote(ssm$model_path),"\n")	
	add_line()

	# model path
	cat("population: ",dQuote(ssm$pop_name),"\n")	
	add_line()

	# state variables
	cat("state variables:",dQuote(ssm$state_variables),"\n")
	add_line()

	# theta
	cat("theta:",dQuote(names(ssm$theta)),"\n")
	add_line()

	# prior
	cat("prior:\n")
	l_ply(ssm$priors, print_prior)
	add_line()

	# theta
	cat("current theta:\n")
	print(ssm$theta %>% signif(digits=3))
	add_line()

	# covmat
	cat("current covariance matrix:\n")
	print(ssm$covmat %>% signif(digits=3))
	add_line()

	# summary
	if(!is.null(ssm$summary)){
		cat("summary:","\n")
		print(ssm$summary)
	}

	# data
	# cat("data:","\n")
	# ssm$data %>% as.tbl %>% print
	# add_line()

	# TODO: reactions, observations, start_date, inputs
}