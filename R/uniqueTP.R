uniqueTIME <- function(object, s, t) {
	return( .Call("uniqueTIME", object, s, t, PACKAGE="TPmsm") )
}

uniqueCOV <- function(object, x) {
	return( .Call("uniqueCOV", object, x, PACKAGE="TPmsm") )
}
