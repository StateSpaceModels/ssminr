#'Convert symbols for latex and html
#'
#'Symbols must be expressed as a string to be converted.
#' @param x character vector containing symbols to be converted in latex or html.
#' @param  to charcater, either \code{"html"} or \code{"latex"}
#' @export
#' @import stringr
convert_symbol <- function(x, to=c("html","latex")) {

	to <- match.arg(to)

	string_greek <- c("alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta", "iota", "kappa", "lambda", "mu", "nu", "xi", "omicron", "pi", "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega", "Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta", "Eta", "Theta", "Iota", "Kappa", "Lambda", "Mu", "Nu", "Xi", "Omicron", "Pi", "Rho", "Sigma", "Tau", "Upsilon", "Phi", "Chi", "Psi", "Omega")

	# make sure it's full name (otherwise issue, e.g. beta -> b+eta)
	replace_symbols <- sprintf("\\b%s\\b", string_greek)

	html_symbols <- sprintf("&%s;", string_greek)
	names(html_symbols) <- replace_symbols

	latex_symbols <- sprintf("\\%s", string_greek)
	names(latex_symbols) <- replace_symbols

	if(to == "html"){
		str_replace_all(x, html_symbols)
	} else if(to == "latex"){
		str_replace_all(x, latex_symbols)
	}

}