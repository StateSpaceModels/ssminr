% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/input.r
\name{input}
\alias{input}
\title{Create SSM input}
\usage{
input(
  name,
  description = NULL,
  value = NULL,
  prior = NULL,
  transformation = NULL,
  to_resource = NULL,
  sde = NULL,
  tag = c("none", "remainder", "pop_size")
)
}
\arguments{
\item{name}{character, name of the input}

\item{description}{character, description of the input}

\item{value}{numeric, define the value of the input. Forced inputs can be specified by a named numeric vector with dates as names.}

\item{prior}{define the prior of the input, as returned by a \code{\link{prior}} helper}

\item{transformation}{define the transformation of the input (see example)}

\item{to_resource}{define the back-transformation of the input. In order to make predictions after fitting your data, specify how to invert the transformation relation at a later time than t0.}

\item{sde}{define a stochastic differential equation on the input (see example)}

\item{tag}{character, tag for specific inputs. Set to one among:
\itemize{
    \item "remainder" if the population size is assumed constant, the tagged state variable will be used as a remainder (see example)
    \item "pop_size"  if the population size is assumed constant, the tagged parameter will be used to set the population size
}}
}
\description{
Create an input parameter for SSM
}
\examples{
\dontrun{
 TODO
}
}
\seealso{
\code{\link{prior}}
}
