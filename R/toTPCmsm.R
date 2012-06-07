toTPCmsm <- function(lst, UT, UX, s, t, x, statenames) {
	return( .Call("toTPCmsm", lst, UT, UX, s, t, x, statenames, PACKAGE="TPmsm") )
}
