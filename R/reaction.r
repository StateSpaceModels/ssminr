#'Create SSM reaction
#'
#'Create a reaction for SSM
#' @param from,to character, name of starting and ending state of the reaction
#' @param description character, description of the reaction
#' @param rate character, rate of the reaction
#' @param accumulators character vector, one or more accumulators. In \code{reaction_split}, several accumulators can be defined for each split reaction using a \code{list} of vectors.
#' @param keywords character vector
#' @export
#' @name reaction
#' @aliases reaction
reaction <- function(from, to, description=NULL, rate, accumulators=NULL, keywords=NULL) {

	reaction <- list(from=from, to=to, description=description, rate=rate)

	if(!is.null(accumulators)){
		reaction$accumulators <- as.list(accumulators)
	}

	if(!is.null(keywords)){
		reaction$keywords <- match.arg(keywords, choices = c("waiting","transmission","linear"), several.ok = TRUE)	
	}

	return(reaction)

}

#' @param split_to a vector whose names will be the \code{to}, and whose values will be the splitting probabilities between the \code{to}.
#' @name reaction
#' @aliases reaction_split
reaction_split <- function(from, split_to, description=NULL, rate, accumulators=NULL, keywords=NULL) {

	if(length(accumulators)>length(split_to)){
		stop("Length of ", sQuote("accumulators"), " and ", sQuote("split_to"), " differ")
	}

	to <- names(split_to)

	reactions <- mapply(reaction, from, to, description, rate, accumulators, MoreArgs = list(keywords=keywords), SIMPLIFY=FALSE, USE.NAMES=FALSE)

	for(i in seq_along(split_to)){

		reactions[[i]]$split <- split_to[[i]]
		# mark reaction order in splitting: if erlang_shape > 1 we need to erlangify only one of the splitted reaction
		reactions[[i]]$keywords <- c(reactions[[i]]$keywords, paste0("split",i))

	}

	return(reactions)

}
