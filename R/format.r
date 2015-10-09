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



#'Simplify symbolic expression
#'
#'Simplify symbolic epxression using rSymPy
#' @param x character vector containing symbolic expressions to be simplified
#' @param  var_names variables appearing in \code{x}
#' @export
#' @import rSymPy
sympy_simplify <- function(x, var_names) {

	# remove space
	x <- str_replace_all(x, " ","")

	# special functions
	reserved <- c('U', 'x', 't', 'E', 'LN2', 'LN10','LOG2E', 'LOG10E', 'PI', 'SQRT1_2', 'SQRT2') #JS Math Global Object
	special_functions <- c('terms_forcing', 'heaviside', 'ramp', 'slowstep', 'sigmoid', 'sin', 'cos', 'correct_rate', 'ssm_correct_rate', 'sqrt', 'pow', 'exp', 'log')

	sink(file.path(tempdir(),"null"))
	all_var <- c(var_names, reserved, special_functions) %>% sapply(Var)
	new_x <- try(sapply(x, function(xx) {print(Var(xx))})) %>% unname
	sink()

	if(inherits(new_x, "try-error")){
		return(x)
	} else {
		return(new_x)
	}

}



