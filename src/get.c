
#include "defines.h"

/*
Author:
	Artur Araujo <artur.stat@gmail.com>

Description:
	This function reads bidimensional data in a way
		identical to how humans read a cartesian plot.
		In this case a step function is considered.

Parameters:
	X[in]			pointer to X first element.
	Y[in]			pointer to Y first element.
	index[in]		pointer to index first element.
	len[in]			pointer to length of X, Y and index.
	i[inout]		pointer to index to start the search.
	a[in]			pointer to abscissa value to read the ordinate at.
	o[out]			pointer to ordinate value.

Return value:
	This function doesn't return a value.

Remarks:
	Vectors X, Y and index must have the same length.
	For this function to work properly vector index
		must indicate the permutation of vector X
		sorted by ascending order.
*/

void getOrdinateI(
	CdoubleCP X,
	CdoubleCP Y,
	CintCP index,
	CintCP len,
	intCP i,
	CdoubleCP a,
	doubleCP o)
{
	if (*i >= *len-1) *o = Y[index[*len-1]];
	else if (X[index[*i]] <= *a) {
		int j = *i;
		*i = (*len-1+j)/2; // the midpoint
		if (X[index[*i+1]] > *a) *i = j;
		for (; *i < *len-1; (*i)++) {
			if (X[index[*i+1]] > *a) break;
		}
		*o = Y[index[*i]];
	}
	return;
} // getOrdinateI

/*
Author:
	Artur Araujo <artur.stat@gmail.com>

Description:
	Computes first e index of a vector T where T[index[e]] > t.

Parameters:
	T[in]			pointer to T first element.
	index[in]		pointer to index first element.
	t[in]			pointer to t value.
	len[in]			pointer to length of T and index.
	i[in]			pointer to index to start the search.
	e[out]			pointer to index found.

Return value:
	This function doesn't return a value.

Remarks:
	Vector index must indicate the permutation
		of vector T sorted by ascending order.
*/

void getIndexI(
	CdoubleCP T,
	CintCP index,
	CdoubleCP t,
	CintCP len,
	CintCP i,
	intCP e)
{
	if (*i >= *len) *e = *len;
	else {
		if (*i < 0) {
			*e = (*len-1)/2; // the midpoint
			if (T[index[*e]] > *t) *e = 0;
		} else {
			*e = (*len-1+*i)/2; // the midpoint
			if (T[index[*e]] > *t) *e = *i;
		}
		for (; *e < *len; (*e)++) {
			if (T[index[*e]] > *t) break; // find index
		}
	}
	return;
} // getIndexI

/*
Author:
	Artur Araujo <artur.stat@gmail.com>

Description:
	Computes last e index of a vector T where T[index[e]] <= t.
	The search is done backwards.

Parameters:
	T[in]			pointer to T first element.
	index[in]		pointer to index first element.
	t[in]			pointer to t value.
	len[in]			pointer to length of T and index.
	i[in]			pointer to index to start the search.
	e[out]			pointer to index found.

Return value:
	This function doesn't return a value.

Remarks:
	Vector index must indicate the permutation
		of vector T sorted by ascending order.
*/

void getBackIndexI(
	CdoubleCP T,
	CintCP index,
	CdoubleCP t,
	CintCP len,
	CintCP i,
	intCP e)
{
	if (*i < 0) *e = -1;
	else {
		if (*i < *len) {
			*e = (*len-1+*i)/2; // the midpoint
			if (T[index[*e]] < *t) *e = *i;
		} else *e = *len-1;
		for (; *e >= 0; (*e)--) {
			if (T[index[*e]] <= *t) break; // find index
		}
	}
	return;
} // getBackIndexI
