#'Create SSM input
#'
#'Create an input parameter for SSM
#' @param name character, name of the input
#' @param description character, description of the input
#' @param value numeric or character, define the value of the input
#' @param prior define the prior of the input, as returned by a \code{\link{prior}} helper
#' @param transformation define the transformation of the input (see example)
#' @param sde define a stochastic differential equation on the input (see example)
#' @param file_name character, name of the file 
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
input <- function(name, description=NULL, value=NULL, prior=NULL, transformation=NULL, to_resource=NULL, sde=NULL, force=FALSE, tag=c("none","remainder","pop_size")) {

	tag <- match.arg(tag)

	require <- NULL

	# TODO: if value is a vector with dates as names => force time varying value
	# need to find a way to write the file outside/inside this function

	# TODO do some check here
	if(!is.null(names(value))){

		for(x in names(value)){

			if(is.null(prior[[x]])){

				if(force){
					stop("Forcing for multiple parameters doesn't work yet")
				} else {
					# define dirac on value
					prior[[x]] <- dirac(value[[x]])
				}
			}

		}

	} else if(!is.null(value) && is.null(prior)){

		if(force){
			require <- list(path=file.path("data",sprintf("%s.csv",name)), fields=c("date",name))				
		} else {

			# define dirac on value
			prior <- dirac(value)
		}
	}

	list(name=name, description=description, value=value, prior=prior, transformation=transformation, to_resource=to_resource, sde=sde, require=require, tag=tag)
	
}