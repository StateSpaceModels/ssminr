remove_null <- function(x) {

	if(is.list(x)){
		x <- x[!sapply(x,is.null)]
	}

	return(x)
}

clean_args <- function(x) {

	x <- x %>% remove_null %>% .[!sapply(.,is.logical) | as.logical(.)] %>% unlist %>% str_replace("TRUE","")
	
	return(x)
}


get_name <- function(x) {

	if(is.list(x)){
		return(sapply(x, function(xx) {xx$name}))
	}

}

# get_args <- function(args_names) {

# 	browser()

# 	arg_names <- formals(fun=sys.function(which=2)) %>% names

# 	env <- parent.frame()

# 	args <- sapply(arg_names, get, envir=env)

# 	return(args)

# }