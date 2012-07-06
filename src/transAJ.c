/***********************************************/
/*** AALEN-JOHANSEN TRANSITION PROBABILITIES ***/
/***********************************************/

#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <Rdefines.h>
#include "defines.h"
#include "boot.h"
#include "get.h"
#include "sort.h"

#define saveAJ \
	P[*b+*nb*i] = p[0]; \
	if (p[2] == 0) P[*b+*nb*(i+*nt)] = 0; \
	else P[*b+*nb*(i+*nt)] = p[2]*p[1]; \
	P[*b+*nb*(i+*nt*2)] = p[2]; \
	i++; \

#define transAJ22I \
	for (; j < eS; j++) { \
		if (S[index1[j]] > UT[i]) { \
			saveAJ \
		} \
		if (T1[index1[j]] < S[index1[j]] && E[index1[j]]) { \
			for (y = j+1, nr = 1; y < *len; y++) nr += (T1[index1[y]] < S[index1[j]]); \
			p[2] *= 1-(double)E[index1[j]]/nr; \
		} \
	} \

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes the transition probabilities:
		p11(s,t) = P(Z>t|Z>s) = P(Z>t)/P(Z>s)
		p12(s,t) = P(Z<=t,T>t|Z>s) = P(s<Z<=t,T>t)/P(Z>s)
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
	P[out]			pointer to a (nb)x(nt)x3 probability array.
	b[in]			pointer to row index.

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
	double P[*nb*(*nt)*3],
	CintCP b)
{
	register int i, y;
	int j = 0, eT, eS, nr;
	double p[3] = {1, 0, 1};
	getIndexI(T1, index0, &UT[0], len, &j, &eT); // determine first index
	int k = eT; // backup first index
	getIndexI(T1, index0, &UT[*nt-1], len, &k, &eT); // determine last index
	getIndexI(S, index1, &UT[0], len, &j, &eS); // determine first index
	for (j = eS, i = 0; k < eT; k++) { // loop through the sample until last index is reached
		getIndexI(S, index1, &T1[index0[k]], len, &j, &eS); // determine index
		transAJ22I // compute transition probability
		if (T1[index0[k]] > UT[i]) {
			saveAJ // save transition probabilities
		}
		p[1] += p[0]*(T1[index0[k]] < S[index0[k]])/( p[2]*(*len-k) ); // compute transition probability
		p[0] *= 1-(double)E1[index0[k]]/(*len-k); // compute transition probability
	}
	getIndexI(S, index1, &UT[*nt-1], len, &j, &eS); // determine index
	transAJ22I // compute transition probability
	saveAJ // save transition probabilities
	for (; i < *nt; i++) { // needed for bootstrap
		P[*b+*nb*i] = P[*b+*nb*(i-1)];
		P[*b+*nb*(i+*nt)] = P[*b+*nb*(i-1+*nt)];
		P[*b+*nb*(i+*nt*2)] = P[*b+*nb*(i-1+*nt*2)];
	}
	return;
} // transAJI

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes a transition probability vector based
		on the Aalen-Johansen estimator.

Parameters:
	object			an object of class 'AJ'.
	UT				unique times vector.
	nboot			number of bootstrap samples.

Return value:
	Returns a list where the first element is a
		(nboot)x(nt)x3 array of transition probabilities,
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
	PROTECT( P = alloc3DArray(REALSXP, *INTEGER(nboot), nt, 3) );
	PROTECT( list = NEW_LIST(2) );
	#ifdef _OPENMP
	#pragma omp parallel if(*INTEGER(nboot) > 1)
	#endif
	{
		int b;
		int *index0 = (int*)malloc( len*sizeof(int) ); // allocate memory block
		int *index1 = (int*)malloc( len*sizeof(int) ); // allocate memory block
		double *WORK0 = (double*)malloc( len*sizeof(double) ); // allocate memory block
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
			transAJI(&len, REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), index0, index1, &nt, REAL(UT), INTEGER(nboot), REAL(P), &b); // compute transition probabilities
		}
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for (b = 1; b < *INTEGER(nboot); b++) {
			boot_ii(index0, index1, &len); // bootstrap indexes
			order_di(REAL(T1), INTEGER(E1), index0, len, FALSE, FALSE, TRUE, WORK0, WORK1); // get permuation
			order_di(REAL(S), INTEGER(E), index1, len, FALSE, FALSE, TRUE, WORK0, WORK1); // get permuation
			transAJI(&len, REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), index0, index1, &nt, REAL(UT), INTEGER(nboot), REAL(P), &b); // compute transition probabilities
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
