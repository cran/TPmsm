
#include <Rdefines.h>
#include "defines.h"
#include "sort.h"

/*
Author:
	Artur Araujo <artur.stat@gmail.com>

Description:
	Concatenates vectors T1 and S and then computes the unique
		times between a given time s and a given time t.

Parameters:
	len[in]			pointer to length of T1 and S.
	T1[in]			pointer to T1 first element.
	S[in]			pointer to S first element.
	s[in]			pointer to s time.
	t[in]			pointer to t time.
	n[inout]		pointer to length of UT.
	UT[out]			pointer to unique times vector.

Return value:
	This function doesn't return a value.
*/

static void uniqueTS(
	CintCP len,
	Cdouble T1[*len],
	Cdouble S[*len],
	CdoubleCP s,
	CdoubleCP t,
	intCP n,
	double UT[*n])
{
	register int i, e;
	for (i = 0; i < *len; i++) { // concatenate vectors T1 and S
		UT[i+1] = T1[i]; // note that UT[0] is reserved for *s
		UT[*len+i+1] = S[i];
	}
	UT[0] = *s;
	sort_d(UT, *n, FALSE, FALSE); // sort vector
	for (i = 0; i < *n; i++) {
		if (UT[i] >= *s) break; // determine first index
	}
	for (e = i; e < *n; e++) {
		if (UT[e] > *t) break; // determine last index
	}
	for (UT[0] = UT[i], i++, *n = 1; i < e; i++) {
		if (UT[i] != UT[i-1]) UT[(*n)++] = UT[i]; // compute unique vector
	}
	return;
} // uniqueTS

/*
Author:
	Artur Araujo <artur.stat@gmail.com>

Description:
	Concatenates X and x and computes it's unique vector.

Parameters:
	len[in]			pointer to length of X.
	X[in]			pointer to X first element.
	nx[in]			pointer to length of x.
	x[in]			pointer to x first element.
	UX[out]			pointer to unique vector.
	n[out]			pointer to length of UX.

Return value:
	This function doesn't return a value.
*/

static void uniqueX(
	CintCP len,
	Cdouble X[*len],
	CintCP nx,
	Cdouble x[*nx],
	double UX[*len+*nx],
	intCP n)
{
	register int i;
	for (i = 0; i < *len; i++) UX[i] = X[i];
	for (i = 0; i < *nx; i++) UX[*len+i] = x[i];
	sort_d(UX, *len+*nx, FALSE, FALSE); // sort vector
	for (i = 1, *n = 1; i < *len+*nx; i++) {
		if (UX[i] != UX[i-1]) UX[(*n)++] = UX[i]; // compute unique vector
	}
	return;
} // uniqueX

/*
Author:
	Artur Araujo <artur.stat@gmail.com>

Description:
	Concatenates vectors T1 and S and then computes the unique
		times between a given time s and a given time t.

Parameters:
	object			an object of class 'survTP' or 'survTPC'.
	s			first time value to compute the probability at.
	t			second time value to compute the probability at.

Return value:
	Returns a vector of unique times.
*/

SEXP uniqueTIME(
	SEXP object,
	SEXP s,
	SEXP t)
{
	SEXP data, T1, S;
	data = VECTOR_ELT(object, 0);
	T1 = VECTOR_ELT(data, 0);
	S = VECTOR_ELT(data, 2);
	int len = GET_LENGTH(T1), n = 2*len+1;
	SEXP UT;
	PROTECT( UT = NEW_NUMERIC(n) );
	uniqueTS( &len, REAL(T1), REAL(S), REAL(s), REAL(t), &n, REAL(UT) ); // compute unique times
	SET_LENGTH(UT, n);
	UNPROTECT(1);
	return UT;
} // uniqueTIME

/*
Author:
	Artur Araujo <artur.stat@gmail.com>

Description:
	Concatenates X and x and computes it's unique vector.

Parameters:
	object			an object of class 'survTPC'.
	x			single covariate value.

Return value:
	Returns a vector of unique covariate values.
*/

SEXP uniqueCOV(
	SEXP object,
	SEXP x)
{
	SEXP data, X;
	data = VECTOR_ELT(object, 0);
	X = VECTOR_ELT(data, 4);
	int len, nx, n;
	len = GET_LENGTH(X);
	nx = GET_LENGTH(x);
	n = len+nx;
	SEXP UX;
	PROTECT( UX = NEW_NUMERIC(n) );
	uniqueX(&len, REAL(X), &nx, REAL(x), REAL(UX), &n); // compute unique times
	SET_LENGTH(UX, n);
	UNPROTECT(1);
	return UX;
} // uniqueCOV
