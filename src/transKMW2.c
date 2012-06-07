/*************************************/
/*** KMW2 TRANSITION PROBABILITIES ***/
/*************************************/

#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <Rdefines.h>
#include "defines.h"
#include "boot.h"
#include "get.h"
#include "sort.h"

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes the transition probabilities:
		p11(s,t) = P(Z>t|Z>s) = P(Z>t)/P(Z>s)
		p13(s,t) = P(T<=t|Z>s) = P(Z>s,T<=t)/P(Z>s)
		p23(s,t) = P(T<=t|Z<=s,T>s) = P(Z<=s,s<T<=t)/P(Z<=s,T>s)

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

static void transKMW3I(
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
	register int i;
	int j = 0, e;
	double aux[4], p[2];
	getIndexI(T1, index0, &UT[0], len, &j, &e); // determine first index
	for (aux[0] = 1, aux[3] = 1; j < e; j++) { // loop through the sample until last index is reached
		aux[2] = (double)E1[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[3] -= aux[2]; // compute probability
	}
	getIndexI(T1, index0, &UT[*nt-1], len, &j, &e); // determine last index
	for (i = 0; j < e; j++) { // loop through the sample until last index is reached
		while (T1[index0[j]] > UT[i]) P[*b+*nb*i++] = aux[3]; // save transition probability
		aux[2] = (double)E1[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[3] -= aux[2]; // compute probability
	}
	for (; i < *nt; i++) P[*b+*nb*i] = aux[3];
	j = 0;
	getIndexI(S, index1, &UT[0], len, &j, &e); // determine first index
	for (aux[0] = 1, aux[3] = 0; j < e; j++) { // loop through the sample until last index is reached
		aux[2] = (double)E[index1[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[3] += aux[2]; // compute probability
	}
	getIndexI(S, index1, &UT[*nt-1], len, &j, &e); // determine last index
	for (i = 0, p[0] = 0, p[1] = 0; j < e; j++) { // loop through the sample until last index is reached
		while (S[index1[j]] > UT[i]) {
			P[*b+*nb*(i+*nt)] = p[0]; // save probability
			P[*b+*nb*(i+*nt*2)] = p[1]; // save probability
			i++;
		}
		aux[2] = (double)E[index1[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		if (T1[index1[j]] <= UT[0]) p[1] += aux[2]; // compute probability
		else p[0] += aux[2]; // compute probability
	}
	for (; i < *nt; i++) {
		P[*b+*nb*(i+*nt)] = p[0];
		P[*b+*nb*(i+*nt*2)] = p[1];
	}
	for (i = *nt-1; i >= 0; i--) {
		P[*b+*nb*(i+*nt)] /= P[*b];
		P[*b+*nb*(i+*nt*2)] /= 1-P[*b]-aux[3];
		P[*b+*nb*i] /= P[*b];
	}
	return;
} // transKMW3I

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes the transition probabilities:
		p11(s,t) = P(Z>t|Z>s) = P(Z>t)/P(Z>s)
		p13(s,t) = P(T<=t|Z>s) = P(Z>s,T<=t)/P(Z>s)
		p23(s,t) = P(T<=t|Z<=s,T>s) = P(Z<=s,s<T<=t)/P(Z<=s,T>s)

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

static void transKMW4I(
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
	register int i;
	int j = 0, e;
	double aux[4], p[2];
	getIndexI(T1, index0, &UT[0], len, &j, &e); // determine first index
	for (aux[0] = 1, aux[3] = 0; j < e; j++) { // loop through the sample until last index is reached
		aux[2] = (double)E1[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[3] += aux[2]; // compute probability
	}
	getIndexI(T1, index0, &UT[*nt-1], len, &j, &e); // determine last index
	for (i = 0; j < e; j++) { // loop through the sample until last index is reached
		while (T1[index0[j]] > UT[i]) P[*b+*nb*i++] = aux[3]; // save transition probability
		aux[2] = (double)E1[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[3] += aux[2]; // compute probability
	}
	for (; i < *nt; i++) P[*b+*nb*i] = aux[3];
	for (; j < *len; j++) { // loop through the sample until last index is reached
		aux[2] = (double)E1[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[3] += aux[2]; // compute probability
	}
	for (i = 0; i < *nt; i++) P[*b+*nb*i] = aux[3]-P[*b+*nb*i];
	j = 0;
	getIndexI(S, index1, &UT[0], len, &j, &e); // determine first index
	for (aux[0] = 1; j < e; j++) { // loop through the sample until last index is reached
		aux[1] = 1-(double)E[index1[j]]/(*len-j); // compute needed factor
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
	}
	getIndexI(S, index1, &UT[*nt-1], len, &j, &e); // determine last index
	for (i = 0, p[0] = 0, p[1] = 0; j < e; j++) { // loop through the sample until last index is reached
		while (S[index1[j]] > UT[i]) {
			P[*b+*nb*(i+*nt)] = p[0]; // save probability
			P[*b+*nb*(i+*nt*2)] = p[1]; // save probability
			i++;
		}
		aux[2] = (double)E[index1[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		if (T1[index1[j]] <= UT[0]) p[1] += aux[2]; // compute probability
		else p[0] += aux[2]; // compute probability
	}
	for (; i < *nt; i++) {
		P[*b+*nb*(i+*nt)] = p[0];
		P[*b+*nb*(i+*nt*2)] = p[1];
	}
	for (; j < *len; j++) { // loop through the sample until last index is reached
		aux[2] = (double)E[index1[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		if (T1[index1[j]] <= UT[0]) p[1] += aux[2]; // compute probability
	}
	for (i = *nt-1; i >= 0; i--) {
		P[*b+*nb*(i+*nt)] /= P[*b];
		P[*b+*nb*i] /= P[*b];
		P[*b+*nb*(i+*nt*2)] /= p[1];
	}
	return;
} // transKMW4I

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes a transition probability vector based
		on the presmoothed Kaplan-Meier weights estimator.

Parameters:
	object			an object of class 'KMW2'.
	UT				unique times vector.
	nboot			number of bootstrap samples.
	methodest		an integer indicating the desired method.

Return value:
	Returns a list where the first element is a
		(nboot)x(nt)x3 array of transition probabilities,
		and the second element is NULL.
*/

SEXP TransPROBKMW2(
	SEXP object,
	SEXP UT,
	SEXP nboot,
	SEXP methodest)
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
	void (*func)(CintCP, CdoubleCP, CintCP, CdoubleCP, CintCP, CintCP, CintCP, CintCP, CdoubleCP, CintCP, doubleCP, CintCP);
	switch ( *INTEGER(methodest) ) {
		case 3:
			func = transKMW3I;
			break;
		default:
			func = transKMW4I;
	}
	#ifdef _OPENMP
	#pragma omp parallel if(*INTEGER(nboot) > 1)
	#endif
	{
		int b;
		int *index0 = (int*)malloc( len*sizeof(int) ); // allocate memory block
		int *index1 = (int*)malloc( len*sizeof(int) ); // allocate memory block
		double *WORK0 = (double*)malloc( len*sizeof(double) ); // allocate memory block
		int *WORK1 = (int*)malloc( len*sizeof(int) ); // allocate memory block
		unsigned int iseed = (unsigned int)time(NULL) ^ (unsigned int)omp_get_thread_num(); // save per thread seed
		srand(iseed); // set seed
		#ifdef _OPENMP
		#pragma omp single
		#endif
		{
			b = 0;
			indx_ii(index0, index1, &len); // initialize indexes
			order_di(REAL(T1), INTEGER(E1), index0, len, FALSE, FALSE, TRUE, WORK0, WORK1); // get permuation
			order_di(REAL(S), INTEGER(E), index1, len, FALSE, FALSE, TRUE, WORK0, WORK1); // get permuation
			func(&len, REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), index0, index1, &nt, REAL(UT), INTEGER(nboot), REAL(P), &b); // compute transition probabilities
		}
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for (b = 1; b < *INTEGER(nboot); b++) {
			boot_ii(index0, index1, &len); // bootstrap indexes
			order_di(REAL(T1), INTEGER(E1), index0, len, FALSE, FALSE, TRUE, WORK0, WORK1); // get permuation
			order_di(REAL(S), INTEGER(E), index1, len, FALSE, FALSE, TRUE, WORK0, WORK1); // get permuation
			func(&len, REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), index0, index1, &nt, REAL(UT), INTEGER(nboot), REAL(P), &b); // compute transition probabilities
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
} // TransPROBKMW2
