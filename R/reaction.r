#'Create SSM reaction
#'
#'Create a reaction for SSM
#' @param from,to character, name of starting and ending state of the reaction
#' @param description character, description of the reaction
#' @param rate character, rate of the reaction
#' @param accumulators character vector, one or more accumulators
#' @export
#' @name reaction
#' @aliases reaction
reaction <- function(from, to, description=NULL, rate, accumulators=NULL) {

	reaction <- list(from=from, to=to, description=description, rate=rate)

	if(!is.null(accumulators)){
		reaction$accumulators <- as.list(accumulators)
	}

	return(reaction)

}

#' @param split_to a vector whose names will be the \code{to}, and whose values will be the splitting probabilities between the \code{to}.
#' @name reaction
#' @aliases reaction_split
reaction_split <- function(from, split_to, description=NULL, rate, accumulators=NULL) {

	
	to <- names(split_to)
	rate <- paste(protect(split_to), protect(rate), sep="*")

	return(mapply(reaction, from, to, description, rate, accumulators, SIMPLIFY=FALSE, USE.NAMES=FALSE))

}