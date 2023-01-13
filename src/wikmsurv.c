
#include "defines.h"

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Computes vector of conditional Kaplan-Meier survival
    probabilities of the censoring variable (C).
    G(t|x) = P(C>t|X=x)

Parameters:
  len[in]           pointer to length of T, E, K, index and SV.
  T[in]             pointer to T first element.
  E[in]             pointer to E first element.
  K[in]             pointer to K first element
  index[in]         pointer to index first element.
  end[in]           pointer to last index of T.
  SV[out]           pointer to conditional survival probabilities vector.

Return value:
  This function doesn't return a value.

Remarks:
  Vector index must indicate the permutation of vector T
    sorted by ascending order.
  If there is no need to have the values of vector K saved
    for later use, it is safe to pass the very same vector
    K as vector SV. This way memory usage is lower.
*/

void wikmsurv(
	CintCP len,
	Cdouble T[*len],
	Cint E[*len],
	Cdouble K[*len],
	Cint index[*len],
	CintCP end,
	double SV[*len])
{
	register int i, j;
	double n, r, d;
	for (i = *len-1, n = 0; i >= *end; i--) { // loop through the sample backwards until last index is reached
		n += K[index[i]]; // weight the individuals in risk
	}
	while (i >= 0) { // loop through the sample backwards until zero index is reached
		n += K[index[i]]; // weight the individuals in risk
		r = E[index[i]]*K[index[i]]; // initialize weighting of the events
		d = (1-E[index[i]])*K[index[i]]; // initialize weighting of the censures
		for (j = i, i--; i >= 0 && T[index[i]] == T[index[i+1]]; i--) { // loop backwards until time changes or zero index is reached
			n += K[index[i]]; // weight the individuals in risk
			r += E[index[i]]*K[index[i]]; // weight the events
			d += (1-E[index[i]])*K[index[i]]; // weight the censures
		}
		for (;j > i+1; j--) SV[index[j]] = 1; // save factor
		if (n-r != 0) SV[index[j]] = 1-d/(n-r); // save factor
		else SV[index[j]] = 1; // save factor
	}
	for (i = 1; i < *end; i++) {
		if (T[index[i]] == T[index[i-1]] && SV[index[i]] != 1) continue; // write only once to an index
		SV[index[i]] *= SV[index[i-1]]; // compute conditional survival probabilities
	}
	return;
} // wikmsurv
