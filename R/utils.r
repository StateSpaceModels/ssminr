remove_null <- function(x) {

	if(is.list(x)){
		x <- x[!sapply(x,is.null)]
	}

	return(x)
}


get_name <- function(x) {

	if(is.list(x)){
		return(sapply(x, function(xx) {xx$name}))
	}

}