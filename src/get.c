
#include "defines.h"

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

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
	if (X[index[*i]] > *a) return;
	else if (*i < *len-1) {
		int j = *i;
		*i = (*len-1+j)/2; // the midpoint
		if (X[index[*i+1]] > *a) *i = j;
		for (; *i < *len-1; (*i)++) if (X[index[*i+1]] > *a) break;
	}
	*o = Y[index[*i]];
	return;
} // getOrdinateI

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes first e index of a vector T where T[index[e]] > t.

Parameters:
	T[in]			pointer to T first element.
	index[in]		pointer to index first element.
	t[in]			pointer to t value.
	len[in]			pointer to length of T and index.
	i[in]			pointer to first index.
	e[out]			pointer to last index.

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
		*e = (*len-1+*i)/2; // the midpoint
		if (T[index[*e]] > *t) *e = *i;
		for (; *e < *len; (*e)++) if (T[index[*e]] > *t) break; // determine last index
	}
	return;
} // getIndexI
