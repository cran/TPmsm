AddClasses <- function(object, value) {
	invisible( .Call("SetClass", object, c(class(object)[1], value), PACKAGE="TPmsm") )
}

RemoveClasses <- function(object) {
	invisible( .Call("SetClass", object, class(object)[1], PACKAGE="TPmsm") )
}
