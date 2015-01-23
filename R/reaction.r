#'Create SSM reaction
#'
#'Create a reaction for SSM
#' @param from,to character, name of starting and ending state of the reaction
#' @param description character, description of the reaction
#' @param rate character, rate of the reaction
#' @param accumulators character vector, one or more accumulators
#' @export
reaction <- function(from, to, description=NULL, rate, accumulators=NULL) {

	reaction <- list(from=from, to=to, description=description, rate=rate)

	if(!is.null(accumulators)){
		reaction$accumulators <- as.list(accumulators)
	}

	return(reaction)

}
