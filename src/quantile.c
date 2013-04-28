
#include <R_ext/Arith.h>
#include <R_ext/Boolean.h>
#include "defines.h"
#include "sort.h"

/*
Author:
	Artur Agostinho Araujo <artur.stat@gmail.com>

Description:
	Computes quantiles of a given vector.

Parameters:
	len[in]			pointer to length of V.
	V[inout]		pointer to V first element.
	n[in]			pointer to length of P and Q.
	P[in]			pointer to P first element.
	Q[out]			pointer to Q first element.

Return value:
	This function doesn't return a value.

Remarks:
	Vector V is sorted inside of this function.
*/

void quantile_d(
	CintCP len,
	double V[*len],
	CintCP n,
	Cdouble P[*n],
	double Q[*n])
{
	register int i;
	int j, k;
	double g;
	sort_d(V, *len, FALSE, FALSE); // sort vector
	for (k = 0; k < *len; k++) if ( !ISNAN(V[k]) ) break; // find first NaN or NA
	for (i = 0; i < *n; i++) {
		g = P[i]*(*len-k-1);
		j = g;
		if (j == *len-k-1) Q[i] = V[*len-1]; // compute quantile
		else {
			g -= j;
			Q[i] = (1-g)*V[j+k]+g*V[j+k+1]; // compute quantile
		}
	}
	return;
} // quantile_d
