#'Create SSM input
#'
#'Create an input parameter for SSM
#' @param name character, name of the input
#' @param description character, description of the input
#' @param value numeric, define the value of the input. Forced inputs can be specified by a named numeric vector with dates as names.
#' @param prior define the prior of the input, as returned by a \code{\link{prior}} helper
#' @param transformation define the transformation of the input (see example)
#' @param to_resource define the back-transformation of the input. In order to make predictions after fitting your data, specify how to invert the transformation relation at a later time than t0.
#' @param sde define a stochastic differential equation on the input (see example)
#' @param tag character, tag for specific inputs. Set to one among:
#' \itemize{
#' 	\item "remainder" if the population size is assumed constant, the tagged state variable will be used as a remainder (see example)
#' 	\item "pop_size"  if the population size is assumed constant, the tagged parameter will be used to set the population size
#' }
#' @export
#' @seealso \code{\link{prior}}
#' @examples \dontrun{
#'  TODO
#'}
input <- function(name, description=NULL, value=NULL, prior=NULL, transformation=NULL, to_resource=NULL, sde=NULL, tag=c("none","remainder","pop_size")) {

	tag <- match.arg(tag)

	forced_input <- NULL
	
	if(!is.null(names(value)) && all(!is.na(dates <- as.Date(names(value), format = '%Y-%m-%d')))){

		if(!is.null(prior)){
			stop("Inputs with prior can't be forced.")
		}

		forced_input <- data_frame(dates, value)
		names(forced_input) <- c("date", name)

		# for(x in names(value)){

		# 	if(is.null(prior[[x]])){

		# 		if(forced_input){
		# 			stop("Forcing for multiple parameters doesn't work yet")
		# 		} else {
		# 			# define dirac on value
		# 			prior[[x]] <- dirac(value[[x]])
		# 		}
		# 	}

		# }

	} else if(length(value) == 1 && !is.null(value) && is.null(prior)){

		# define dirac on value
		prior <- dirac(value)
		
	} else if (length(value) > 1) {

		stop("value must be either an atomic value or a numeric vector with dates as name")

	}

	list(name=name, description=description, value=value, prior=prior, transformation=transformation, to_resource=to_resource, sde=sde, forced_input = forced_input, tag=tag)
	
}