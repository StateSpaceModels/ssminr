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
#' @examples 
#' # define a diffusion on log(beta) with volatility parameter "vol"	
#' input(name="beta", description="effective contact rate",transformation="R0/d_infectious", sde=diffusion(volatility="vol", transformation="log"))
diffusion <- function(volatility, transformation=c("none","log")){

	transformation <- match.arg(transformation)

	list(volatility=volatility, transformation=transformation)

}