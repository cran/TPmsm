
#include <Rdefines.h>

SEXP SetClass(
	SEXP object,
	SEXP value)
{
	SET_CLASS(object, value);
	return R_NilValue;
} // SetClass
