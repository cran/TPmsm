/***********************************************/
/*** AALEN-JOHANSEN TRANSITION PROBABILITIES ***/
/***********************************************/

#include <R_ext/Arith.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdlib.h>
#include <time.h>
#include <Rdefines.h>
#include "defines.h"
#include "boot.h"
#include "get.h"
#include "rthreads.h"
#include "sort.h"

#define saveAJ \
	P[*b+*nb*i] = p[0]; /* p11(s,t) */ \
	P[*b+*nb*(i+*nt)] = p[2]*p[1]; /* p12(s,t) */ \
	P[*b+*nb*(i+*nt*2)] = 1-p[0]-p[2]*p[1]; /* p13(s,t) */ \
	if (P[*b+*nb*(i+*nt*2)] < 0) { \
		P[*b+*nb*(i+*nt)] = 1-p[0]; /* p12(s,t) */ \
		P[*b+*nb*(i+*nt*2)] = 0; /* p13(s,t) */ \
	} \
	P[*b+*nb*(i+*nt*3)] = p[2]; /* p22(s,t) */ \
	i++; \

#define transAJ22I \
	for (; y < e[1]; y++) { \
		if (S[index1[y]] > UT[i]) { \
			saveAJ \
		} \
		if (T1[index1[y]] < S[index1[y]] && E[index1[y]]) { \
			for (j = y+1, z = 1; j < *len; j++) z += (T1[index1[j]] < S[index1[y]]); \
			WORK[*len+y] = 1-1.0/z; \
			p[2] *= WORK[*len+y]; \
		} else WORK[*len+y] = 1; \
	} \

/*
Author:
	Artur Agostinho Araujo <artur.stat@gmail.com>

Description:
	Computes the transition probabilities:
		p11(s,t) = P(Z>t|Z>s) = P(Z>t)/P(Z>s)
		p12(s,t) = P(Z<=t,T>t|Z>s) = P(s<Z<=t,T>t)/P(Z>s)
		p13(s,t) = 1-p11(s,t)-p12(s,t)
		p22(s,t) = P(Z<=t,T>t|Z<=s,T>s) = P(Z<=s,T>t)/P(Z<=s,T>s)

Parameters:
	len[in]			pointer to length of T1, E1, S and E.
	T1[in]			pointer to T1 first element.
	E1[in]			pointer to E1 first element.
	S[in]			pointer to S first element.
	E[in]			pointer to E first element.
	index0[in]		pointer to index0 first element.
	index1[in]		pointer to index1 first element.
	nt[in]			pointer to length of UT and number of rows of P.
	UT[in]			pointer to unique times vector.
	nb[in]			pointer to number of rows of P.
	P[out]			pointer to a (nb)x(nt)x4 probability array.
	b[in]			pointer to row index.
	WORK[out]		pointer to WORK first element

Return value:
	This function doesn't return a value.

Remarks:
	Vector index0 must indicate the permutation of vector T1
		sorted by ascending order.
	Vector index1 must indicate the permutation of vector S
		sorted by ascending order.
	Vectors T1, E1, S and E must have the same length.
*/

static void transAJI(
	CintCP len,
	Cdouble T1[*len],
	Cint E1[*len],
	Cdouble S[*len],
	Cint E[*len],
	Cint index0[*len],
	Cint index1[*len],
	CintCP nt,
	Cdouble UT[*nt],
	CintCP nb,
	double P[*nb*(*nt)*4],
	CintCP b,
	double WORK[*len*2])
{
	register int i, j;
	int k, x, y, z, s[2], e[2] = {0, 0};
	double p[3] = {1, 0, 1};
	getIndexI(T1, index0, &UT[0], len, &e[0], &s[0]); // determine first index
	getIndexI(T1, index0, &UT[*nt-1], len, &s[0], &e[0]); // determine last index
	getIndexI(S, index1, &UT[0], len, &e[1], &s[1]); // determine first index
	for (x = s[0], y = s[1], i = 0; x < e[0]; x++) { // loop through the sample until last index is reached
		getIndexI(S, index1, &T1[index0[x]], len, &y, &e[1]); // determine index
		transAJ22I // compute transition probability
		if (T1[index0[x]] > UT[i]) {
			saveAJ // save transition probabilities
		}
		if (T1[index0[x]] < S[index0[x]]) {
			WORK[x] = p[0]/(*len-x); // save transition probability factor
		} else WORK[x] = 0; // save transition probability factor
		p[1] += WORK[x]/p[2]; // compute transition probability
		if (E1[index0[x]]) p[0] *= 1-1.0/(*len-x); // compute transition probability
	}
	getIndexI(S, index1, &UT[*nt-1], len, &y, &e[1]); // determine index
	transAJ22I // compute transition probability
	saveAJ // save transition probabilities
	for (; i < *nt; i++) { // needed for bootstrap
		for (j = 0; j < 4; j++) P[*b+*nb*(i+*nt*j)] = P[*b+*nb*(i-1+*nt*j)];
	}
	for (i--; i >= 0; i--) if ( !ISNAN(P[*b+*nb*(i+*nt)]) ) break; // loop backwards while ISNAN
	for (i++, x = s[0]; i < *nt; i++) {
		getIndexI(T1, index0, &UT[i], len, &x, &e[0]); // determine last index
		x = e[0]; // save index for next search
		if (x > s[0]) break; // break the loop
		P[*b+*nb*(i+*nt)] = 0; // save p12(s,t)
		P[*b+*nb*(i+*nt*2)] = 1-P[*b+*nb*i]; // save p13(s,t)
	}
	for (y = s[1]; i < *nt; i++) {
		getIndexI(T1, index0, &UT[i], len, &x, &e[0]); // determine last index
		x = e[0]; // save index for next search
		getIndexI(S, index1, &UT[i], len, &y, &e[1]); // determine last index
		y = e[1]; // save index for next search
		for (p[1] = 0, p[2] = 1, j = e[0]-1, k = e[1]-1; j >= s[0]; j--) { // loop backwards through the sample until first index is reached
			if (WORK[j] == 0) continue; // don't waste time doing uneeded computations
			getBackIndexI(S, index1, &T1[index0[j]], len, &k, &z); // determine first index
			for (;k > z; k--) p[2] *= WORK[*len+k]; // compute transition probability
			p[1] += WORK[j]*p[2]; // compute transition probability
		}
		P[*b+*nb*(i+*nt)] = p[1]; // save p12(s,t)
		P[*b+*nb*(i+*nt*2)] = 1-P[*b+*nb*i]-p[1]; // save p13(s,t)
		if (P[*b+*nb*(i+*nt*2)] < 0) {
			P[*b+*nb*(i+*nt)] = 1-P[*b+*nb*i]; // save p12(s,t)
			P[*b+*nb*(i+*nt*2)] = 0; // save p13(s,t)
		}
	}
	return;
} // transAJI

/*
Author:
	Artur Agostinho Araujo <artur.stat@gmail.com>

Description:
	Computes a transition probability vector based
		on the Aalen-Johansen estimator.

Parameters:
	object			an object of class 'AJ'.
	UT				unique times vector.
	nboot			number of bootstrap samples.

Return value:
	Returns a list where the first element is a
		(nboot)x(nt)x4 array of transition probabilities,
		and the second element is NULL.
*/

SEXP TransPROBAJ(
	SEXP object,
	SEXP UT,
	SEXP nboot)
{
	SEXP data, T1, E1, S, E;
	data = VECTOR_ELT(object, 0);
	T1 = VECTOR_ELT(data, 0);
	E1 = VECTOR_ELT(data, 1);
	S = VECTOR_ELT(data, 2);
	E = VECTOR_ELT(data, 3);
	int len = GET_LENGTH(T1), nt = GET_LENGTH(UT);
	SEXP P, list;
	PROTECT( P = alloc3DArray(REALSXP, *INTEGER(nboot), nt, 4) );
	PROTECT( list = NEW_LIST(2) );
	#ifdef _OPENMP
	#pragma omp parallel if(*INTEGER(nboot) > 1) num_threads(global_num_threads)
	#endif
	{
		int b;
		int *index0 = (int*)malloc( len*sizeof(int) ); // allocate memory block
		int *index1 = (int*)malloc( len*sizeof(int) ); // allocate memory block
		double *WORK0 = (double*)malloc( len*2*sizeof(double) ); // allocate memory block
		int *WORK1 = (int*)malloc( len*sizeof(int) ); // allocate memory block
		#ifdef _OPENMP
		unsigned int iseed = (unsigned int)time(NULL) ^ (unsigned int)omp_get_thread_num(); // save per thread seed
		#else
		unsigned int iseed = (unsigned int)time(NULL);
		#endif
		srand(iseed); // set seed
		#ifdef _OPENMP
		#pragma omp single
		#endif
		{
			b = 0;
			indx_ii(index0, index1, &len); // initialize indexes
			order_di(REAL(T1), INTEGER(E1), index0, len, FALSE, FALSE, TRUE, WORK0, WORK1); // get permuation
			order_di(REAL(S), INTEGER(E), index1, len, FALSE, FALSE, TRUE, WORK0, WORK1); // get permuation
			transAJI(&len, REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), index0, index1, &nt, REAL(UT), INTEGER(nboot), REAL(P), &b, WORK0); // compute transition probabilities
		}
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for (b = 1; b < *INTEGER(nboot); b++) {
			boot_ii(index0, index1, &len); // bootstrap indexes
			order_di(REAL(T1), INTEGER(E1), index0, len, FALSE, FALSE, TRUE, WORK0, WORK1); // get permuation
			order_di(REAL(S), INTEGER(E), index1, len, FALSE, FALSE, TRUE, WORK0, WORK1); // get permuation
			transAJI(&len, REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), index0, index1, &nt, REAL(UT), INTEGER(nboot), REAL(P), &b, WORK0); // compute transition probabilities
		}
		free(index0); // free memory block
		free(index1); // free memory block
		free(WORK0); // free memory block
		free(WORK1); // free memory block
	}
	SET_ELEMENT(list, 0, P);
	SET_ELEMENT(list, 1, R_NilValue);
	UNPROTECT(2);
	return list;
} // TransPROBAJ
