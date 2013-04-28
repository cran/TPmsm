
#include "defines.h"

/*
Author:
	Artur Agostinho Araujo <artur.stat@gmail.com>

Description:
	Computes unique index vector.

Parameters:
	T[in]			pointer to T first element.
	index[in]		pointer to index first element.
	len[in]			pointer to length of index.
	unique[out]		pointer to unique vector.
	u[out]			pointer to length of unique vector.

Return value:
	This function doesn't return a value.

Remarks:
	Vector index must indicate the permutation of vector T
		sorted by ascending order.
*/

void uniqueI(
	CdoubleCP T,
	CintCP index,
	CintCP len,
	intCP unique,
	intCP u)
{
	register int i;
	for (unique[0] = index[0], i = 1, *u = 1; i < *len; i++) {
		if (T[index[i]] != T[index[i-1]]) unique[(*u)++] = index[i];
	}
	return;
} // uniqueI
