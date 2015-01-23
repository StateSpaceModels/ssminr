#' Define SDE on SSM inputs
#'
#'In SSM one can define stochastic differential equations on input parameters. 
#' @param volatility character, name of the volatility parameter
#' @param transformation character, define a transformation for the input parameter (none by default). The diffusion is then defined for the transformed parameter:
#' \itemize{
#' 		\item \code{"log"} log-transform for positive parameters.
#' }
#' @export
#' @name sde
#' @aliases diffusion
#' @example inst/examples/sde-example.r
diffusion <- function(volatility, transformation=c("none","log")){

	transformation <- match.arg(transformation)

	list(volatility=volatility, transformation=transformation)

}