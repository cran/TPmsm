toTPmsm.default <- function(lst, UT, s, t, statenames) {
	return( .Call("toTPmsm1222", lst, UT, s, t, statenames, PACKAGE="TPmsm") )
}

toTPmsm.KMW2 <- function(lst, UT, s, t, statenames) {
	return( .Call("toTPmsm1323", lst, UT, s, t, statenames, PACKAGE="TPmsm") )
}

toTPmsm.KMPW2 <- function(lst, UT, s, t, statenames) {
	return( .Call("toTPmsm1323", lst, UT, s, t, statenames, PACKAGE="TPmsm") )
}
