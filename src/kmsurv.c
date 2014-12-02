
#include "defines.h"

/*
Author:
	Artur Araujo <artur.stat@gmail.com>

Description:
	Computes vector of Kaplan-Meier survival probabilities.
		S(t) = P(T>t)

Parameters:
	len[in]			pointer to length of T, E, index and SV.
	T[in]			pointer to T first element.
	E[in]			pointer to E first element.
	index[in]		pointer to index first element.
	end[in]			pointer to last index of T.
	SV[out]			pointer to survival probabilities vector.

Return value:
	This function doesn't return a value.

Remarks:
	Vector index must indicate the permutation of vector T
		sorted by ascending order.
*/

void kmsurv(
	CintCP len,
	Cdouble T[*len],
	Cint E[*len],
	Cint index[*len],
	CintCP end,
	double SV[*len])
{
	register int i = 0, j;
	int n, d;
	double p = 1;
	while (i < *end) { // loop through the sample until last index is reached
		n = *len-i; // count the individuals in risk
		d = E[index[i]]; // initialize event count
		for (i++; i < *end && T[index[i]] == T[index[i-1]]; i++) { // loop until time changes or last index is reached
			d += E[index[i]]; // count the events
		}
		p *= 1-(double)d/n; // compute survival probability
		for (j = *len-n; j < i; j++) SV[index[j]] = p; // save survival probability
	}
	return;
} // kmsurv
