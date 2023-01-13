/*************************************/
/*** LIN1 TRANSITION PROBABILITIES ***/
/*************************************/

#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdlib.h>
#include <Rdefines.h>
#include "defines.h"
#include "get.h"
#include "RngStream.h"
#include "RngArray.h"
#include "RngBoot.h"
#include "rthreads.h"
#include "sort.h"

#define isurv_body \
	while (j < e) { \
		n = *len-j; \
		r = E1[index0[j]]; \
		d = 1-E1[index0[j]]; \
		for (j++; j < e && T1[index0[j]] == T1[index0[j-1]]; j++) { \
			r += E1[index0[j]]; \
			d += 1-E1[index0[j]]; \
		} \
		if (n-r != 0) surv *= 1-(double)d/(n-r); \
	} \

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Computes the transition probabilities:
    p11(s,t) = P(Z>t|Z>s) = P(Z>t)/P(Z>s)
    p12(s,t) = P(Z<=t,T>t|Z>s) = P(s<Z<=t,T>t)/P(Z>s)
    p13(s,t) = 1-p11(s,t)-p12(s,t)
    p22(s,t) = P(Z<=t,T>t|Z<=s,T>s) = P(Z<=s,T>t)/P(Z<=s,T>s)

Parameters:
  len[in]           pointer to length of T1, E1, S and E.
  T1[in]            pointer to T1 first element.
  E1[in]            pointer to E1 first element.
  S[in]             pointer to S first element.
  E[in]             pointer to E first element.
  index0[in]        pointer to index0 first element.
  index1[in]        pointer to index1 first element.
  nt[in]            pointer to length of UT and number of rows of P.
  UT[in]            pointer to unique times vector.
  nb[in]            pointer to number of rows of P.
  P[out]            pointer to a (nb)x(nt)x4 probability array.
  b[in]             pointer to row index.

Return value:
  This function doesn't return a value.

Remarks:
  Vector index0 must indicate the permutation of vector T1
    sorted by ascending order.
  Vector index1 must indicate the permutation of vector S
    sorted by ascending order.
  Vectors T1, E1, S and E must have the same length.
*/

static void transLIN1I(
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
	CintCP b)
{
	register int i = 0;
	int j = 0, e, n, r, d, sum[2];
	double surv = 1, p[2];
	getIndexI(T1, index0, &UT[0], len, &j, &e); // determine first index
	isurv_body
	sum[0] = e;
	getIndexI(T1, index0, &UT[*nt-1], len, &j, &e); // determine index
	while (j < e) { // loop through the sample until last index is reached
		n = *len-j; // count individuals at risk
		p[0] = n/surv;
		p[1] = j-sum[0];
		while (T1[index0[j]] > UT[i]) {
			P[*b+*nb*i] = p[0]; // save p11(s,t) factor
			P[*b+*nb*(i+*nt)] = p[1]; // save p12(s,t) term
			i++;
		}
		r = E1[index0[j]]; // initialize event count
		d = 1-E1[index0[j]]; // initialize censure count
		for (j++; j < e && T1[index0[j]] == T1[index0[j-1]]; j++) { // loop until time changes or last index is reached
			r += E1[index0[j]]; // count the events
			d += 1-E1[index0[j]]; // count the censures
		}
		if (n-r != 0) surv *= 1-(double)d/(n-r); // compute survival probability
	}
	p[0] = (*len-j)/surv;
	p[1] = j-sum[0];
	for (; i < *nt; i++) { // needed for bootstrap
		P[*b+*nb*i] = p[0]; // save p11(s,t) factor
		P[*b+*nb*(i+*nt)] = p[1]; // save p12(s,t) term
	}
	j = 0;
	getIndexI(S, index1, &UT[0], len, &j, &e); // determine first index
	surv = 1;
	#define E1 E
	#define index0 index1
	#define T1 S
	isurv_body
	#undef E1
	#undef index0
	#undef T1
	sum[0] -= e;
	sum[1] = 0;
	getIndexI(S, index1, &UT[*nt-1], len, &j, &e); // determine index
	i = 0;
	while (j < e) { // loop through the sample until last index is reached
		p[0] = sum[0]/surv;
		while (S[index1[j]] > UT[i]) {
			P[*b+*nb*(i+*nt)] -= sum[1];
			P[*b+*nb*(i+*nt)] /= surv; // compute and save p12(s,t) factor
			P[*b+*nb*(i+*nt*3)] = p[0]; // save p22(s,t) factor
			i++;
		}
		n = *len-j; // count individuals at risk
		r = E[index1[j]]; // initialize event count
		d = 1-E[index1[j]]; // initialize censure count
		if (T1[index1[j]] <= UT[0]) sum[0] -= 1; // compute sum
		else sum[1] += 1; // compute sum
		for (j++; j < e && S[index1[j]] == S[index1[j-1]]; j++) { // loop until time changes or last index is reached
			r += E[index1[j]]; // count the events
			d += 1-E[index1[j]]; // count the censures
			if (T1[index1[j]] <= UT[0]) sum[0] -= 1; // compute sum
			else sum[1] += 1; // compute sum
		}
		if (n-r != 0) surv *= 1-(double)d/(n-r); // compute survival probability
	}
	p[0] = sum[0]/surv;
	for (; i < *nt; i++) { // needed for bootstrap
		P[*b+*nb*(i+*nt)] -= sum[1];
		P[*b+*nb*(i+*nt)] /= surv; // compute and save p12(s,t) factor
		P[*b+*nb*(i+*nt*3)] = p[0]; // save p22(s,t) factor
	}
	for (i = *nt-1; i >= 0; i--) {
		P[*b+*nb*(i+*nt)] /= P[*b]; // compute and save p12(s,t)
		if (P[*b+*nb*(i+*nt)] < 0) P[*b+*nb*(i+*nt)] = 0;
		P[*b+*nb*i] /= P[*b]; // compute and save p11(s,t)
		if (P[*b+*nb*i] > 1) P[*b+*nb*i] = 1;
		P[*b+*nb*(i+*nt*2)] = 1-P[*b+*nb*i]-P[*b+*nb*(i+*nt)]; // compute and save p13(s,t)
		if (P[*b+*nb*(i+*nt*2)] < 0) {
			P[*b+*nb*(i+*nt)] = 1-P[*b+*nb*i]; // compute and save p12(s,t)
			P[*b+*nb*(i+*nt*2)] = 0; // save p13(s,t)
		}
		P[*b+*nb*(i+*nt*3)] /= P[*b+*nb*(*nt*3)]; // compute and save p22(s,t)
		if (P[*b+*nb*(i+*nt*3)] > 1) P[*b+*nb*(i+*nt*3)] = 1;
	}
	return;
} // transLIN1I

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Computes a transition probability vector based
    on the LIN estimator.

Parameters:
  object            an object of class 'LIN1'.
  UT                unique times vector.
  nboot             number of bootstrap samples.

Return value:
  Returns a list where the first element is a
    (nboot)x(nt)x4 array of transition probabilities,
    and the second element is NULL.
*/

SEXP TransPROBLIN1(
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
	int len = GET_LENGTH(T1), nt = GET_LENGTH(UT), t, nth = 1;
	SEXP P, list;
	PROTECT( P = alloc3DArray(REALSXP, *INTEGER(nboot), nt, 4) );
	PROTECT( list = NEW_LIST(2) );
	if (*INTEGER(nboot) > 1) nth = global_num_threads;
	int **index0 = (int**)malloc( (unsigned int)nth*sizeof(int*) ); // allocate memory block
	if (index0 == NULL) error("TransPROBLIN1: No more memory\n");
	int **index1 = (int**)malloc( (unsigned int)nth*sizeof(int*) ); // allocate memory block
	if (index1 == NULL) error("TransPROBLIN1: No more memory\n");
	double **WORK = (double**)malloc( (unsigned int)nth*sizeof(double*) ); // allocate memory block
	if (WORK == NULL) error("TransPROBLIN1: No more memory\n");
	for (t = 0; t < nth; t++) { // allocate per thread memory
		if ( ( index0[t] = (int*)malloc( (unsigned int)len*sizeof(int) ) ) == NULL ) error("TransPROBLIN1: No more memory\n");
		if ( ( index1[t] = (int*)malloc( (unsigned int)len*sizeof(int) ) ) == NULL ) error("TransPROBLIN1: No more memory\n");
		if ( ( WORK[t] = (double*)malloc( (unsigned int)len*sizeof(double) ) ) == NULL ) error("TransPROBLIN1: No more memory\n");
	}
	#ifdef _OPENMP
	#pragma omp parallel num_threads(nth) private(t)
	#endif
	{
		int b;
		#ifdef _OPENMP
		t = omp_get_thread_num();
		#else
		t = 0;
		#endif
		#ifdef _OPENMP
		#pragma omp single
		#endif
		{
			b = 0;
			indx_ii(&len, index0[t], index1[t]); // initialize indexes
			order_d(REAL(T1), index0[t], len, FALSE, FALSE, WORK[t]); // get permuation
			order_d(REAL(S), index1[t], len, FALSE, FALSE, WORK[t]); // get permuation
			transLIN1I(&len, REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), index0[t], index1[t], &nt, REAL(UT), INTEGER(nboot), REAL(P), &b); // compute transition probabilities
		}
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for (b = 1; b < *INTEGER(nboot); b++) {
			boot_ii(RngArray[t], &len, index0[t], index1[t]); // bootstrap indexes
			order_d(REAL(T1), index0[t], len, FALSE, FALSE, WORK[t]); // get permuation
			order_d(REAL(S), index1[t], len, FALSE, FALSE, WORK[t]); // get permuation
			transLIN1I(&len, REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), index0[t], index1[t], &nt, REAL(UT), INTEGER(nboot), REAL(P), &b); // compute transition probabilities
		}
	}
	for (t = nth-1; t >= 0; t--) {
		free(index0[t]); // free memory block
		free(index1[t]); // free memory block
		free(WORK[t]); // free memory block
	}
	free(index0); // free memory block
	free(index1); // free memory block
	free(WORK); // free memory block
	SET_ELEMENT(list, 0, P);
	SET_ELEMENT(list, 1, R_NilValue);
	UNPROTECT(2);
	return list;
} // TransPROBLIN1
