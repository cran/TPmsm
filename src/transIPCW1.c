/**************************************/
/*** IPCW1 TRANSITION PROBABILITIES ***/
/**************************************/

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

#define invsum_body \
	while (j < e) { \
		n = *len-j; \
		r = E1[index0[j]]; \
		d = 1-E1[index0[j]]; \
		for (j++; j < e && T1[index0[j]] == T1[index0[j-1]]; j++) { \
			r += E1[index0[j]]; \
			d += 1-E1[index0[j]]; \
		} \
		if (n-r != 0) surv *= 1-(double)d/(n-r); \
		if (surv > 0) { \
			for (k = *len-n; k < j; k++) sum += E1[index0[k]]/surv; \
		} \
	} \

/*
Author:
	Artur Agostinho Araujo <artur.stat@gmail.com>

Description:
	Computes the transition probabilities:
		p11(s,t) = P(Z>t|Z>s) = {1-P(Z<=t)}/{1-P(Z<=s)}
		p12(s,t) = P(Z<=t,T>t|Z>s) = {P(Z<=t)-P(Z<=s)-P(s<Z<=t,T<=t)}/{1-P(Z<=s)}
		p13(s,t) = 1-p11(s,t)-p12(s,t)
		p22(s,t) = P(Z<=t,T>t|Z<=s,T>s) = {P(Z<=s)-P(Z<=s,T<=t)}/{P(Z<=s)-P(T<=s)}

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

Return value:
	This function doesn't return a value.

Remarks:
	Vector index0 must indicate the permutation of vector T1
		sorted by ascending order.
	Vector index1 must indicate the permutation of vector S
		sorted by ascending order.
	Vectors T1, E1, S and E must have the same length.
*/

static void transIPCW1I(
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
	register int i = 0, k;
	int j = 0, e, n, r, d;
	double sum[4] = {0, 0, 0, 0}, surv = 1, p[2];
	getIndexI(T1, index0, &UT[0], len, &j, &e); // determine first index
	#define sum sum[0]
	invsum_body
	#undef sum
	getIndexI(T1, index0, &UT[*nt-1], len, &j, &e); // determine index
	while (j < e) { // loop through the sample until last index is reached
		p[0] = (*len-sum[1]-sum[0])/(*len-sum[0]); // compute p11(s,t)
		p[1] = sum[1]/(*len-sum[0]); // compute p12(s,t) term
		while (T1[index0[j]] > UT[i]) {
			if (p[0] < 0) P[*b+*nb*i] = 0;
			else P[*b+*nb*i] = p[0]; // save p11(s,t)
			P[*b+*nb*(i+*nt)] = p[1]; // save p12(s,t) term
			i++;
		}
		n = *len-j; // count individuals at risk
		r = E1[index0[j]]; // initialize event count
		d = 1-E1[index0[j]]; // initialize censure count
		for (j++; j < e && T1[index0[j]] == T1[index0[j-1]]; j++) { // loop until time changes or last index is reached
			r += E1[index0[j]]; // count the events
			d += 1-E1[index0[j]]; // count the censures
		}
		if (n-r != 0) surv *= 1-(double)d/(n-r); // compute survival probability
		if (surv > 0) for (k = *len-n; k < j; k++) sum[1] += E1[index0[k]]/surv; // compute sum
	}
	p[0] = (*len-sum[1]-sum[0])/(*len-sum[0]); // compute p11(s,t)
	p[1] = sum[1]/(*len-sum[0]); // compute p12(s,t) term
	for (; i < *nt; i++) { // needed for bootstrap
		if (p[0] < 0) P[*b+*nb*i] = 0;
		else P[*b+*nb*i] = p[0]; // save p11(s,t)
		P[*b+*nb*(i+*nt)] = p[1]; // save p12(s,t) term
	}
	sum[1] = 0;
	j = 0;
	surv = 1;
	getIndexI(S, index1, &UT[0], len, &j, &e); // determine first index
	#define E1 E
	#define index0 index1
	#define T1 S
	#define sum sum[1]
	invsum_body
	#undef E1
	#undef index0
	#undef T1
	#undef sum
	getIndexI(S, index1, &UT[*nt-1], len, &j, &e); // determine index
	i = 0;
	while (j < e) { // loop through the sample until last index is reached
		p[0] = sum[2]/(*len-sum[0]); // compute p12(s,t) term
		p[1] = 1-sum[3]/(sum[0]-sum[1]); // compute p22(s,t)
		while (S[index1[j]] > UT[i]) {
			P[*b+*nb*(i+*nt)] -= p[0]; // compute and save p12(s,t)
			if (P[*b+*nb*(i+*nt)] < 0) P[*b+*nb*(i+*nt)] = 0;
			P[*b+*nb*(i+*nt*2)] = 1-P[*b+*nb*i]-P[*b+*nb*(i+*nt)]; // compute and save p13(s,t)
			if (P[*b+*nb*(i+*nt*2)] < 0) {
				P[*b+*nb*(i+*nt)] = 1-P[*b+*nb*i]; // compute and save p12(s,t)
				P[*b+*nb*(i+*nt*2)] = 0; // save p13(s,t)
			}
			if (p[1] < 0) P[*b+*nb*(i+*nt*3)] = 0;
			else P[*b+*nb*(i+*nt*3)] = p[1]; // save p22(s,t)
			i++;
		}
		n = *len-j; // count individuals at risk
		r = E[index1[j]]; // initialize event count
		d = 1-E[index1[j]]; // initialize censure count
		for (j++; j < e && S[index1[j]] == S[index1[j-1]]; j++) { // loop until time changes or last index is reached
			r += E[index1[j]]; // count the events
			d += 1-E[index1[j]]; // count the censures
		}
		if (n-r != 0) surv *= 1-(double)d/(n-r); // compute survival probability
		if (surv > 0) {
			for (k = *len-n; k < j; k++) {
				if (T1[index1[k]] <= UT[0]) sum[3] += E[index1[k]]/surv; // compute sum
				else sum[2] += E[index1[k]]/surv; // compute sum
			}
		}
	}
	p[0] = sum[2]/(*len-sum[0]); // compute p12(s,t) term
	p[1] = 1-sum[3]/(sum[0]-sum[1]); // compute p22(s,t)
	for (; i < *nt; i++) { // needed for bootstrap
		P[*b+*nb*(i+*nt)] -= p[0]; // compute and save p12(s,t)
		if (P[*b+*nb*(i+*nt)] < 0) P[*b+*nb*(i+*nt)] = 0;
		P[*b+*nb*(i+*nt*2)] = 1-P[*b+*nb*i]-P[*b+*nb*(i+*nt)]; // compute and save p13(s,t)
		if (P[*b+*nb*(i+*nt*2)] < 0) {
			P[*b+*nb*(i+*nt)] = 1-P[*b+*nb*i]; // compute and save p12(s,t)
			P[*b+*nb*(i+*nt*2)] = 0; // save p13(s,t)
		}
		if (p[1] < 0) P[*b+*nb*(i+*nt*3)] = 0;
		else P[*b+*nb*(i+*nt*3)] = p[1]; // save p22(s,t)
	}
	return;
} // transIPCW1I

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

Return value:
	This function doesn't return a value.

Remarks:
	Vector index0 must indicate the permutation of vector T1
		sorted by ascending order.
	Vector index1 must indicate the permutation of vector S
		sorted by ascending order.
	Vectors T1, E1, S and E must have the same length.
*/

static void transIPCW2I(
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
	register int i = 0, k, y = 0;
	int j = 0, e, n, r, d;
	double sum = 0, surv = 1;
	getIndexI(T1, index0, &UT[0], len, &j, &e); // determine first index
	invsum_body
	getIndexI(T1, index0, &UT[*nt-1], len, &j, &e); // determine index
	while (j < e) { // loop through the sample until last index is reached
		while (T1[index0[j]] > UT[i]) {
			P[*b+*nb*i] = sum; // save p11(s,t) term
			P[*b+*nb*(i+*nt)] = 0; // initialize p12(s,t)
			P[*b+*nb*(i+*nt*3)] = 0; // initialize p22(s,t)
			i++;
		}
		n = *len-j; // count individuals at risk
		r = E1[index0[j]]; // initialize event count
		d = 1-E1[index0[j]]; // initialize censure count
		for (j++; j < e && T1[index0[j]] == T1[index0[j-1]]; j++) { // loop until time changes or last index is reached
			r += E1[index0[j]]; // count the events
			d += 1-E1[index0[j]]; // count the censures
		}
		if (n-r != 0) surv *= 1-(double)d/(n-r); // compute survival probability
		if (surv > 0) for (k = *len-n; k < j; k++) sum += E1[index0[k]]/surv; // compute sum
	}
	for (; i < *nt; i++) { // needed for bootstrap
		P[*b+*nb*i] = sum; // save p11(s,t) term
		P[*b+*nb*(i+*nt)] = 0; // initialize p12(s,t)
		P[*b+*nb*(i+*nt*3)] = 0; // initialize p22(s,t)
	}
	#define e *len
	invsum_body
	#undef e
	for (i = 0; i < *nt; i++) P[*b+*nb*i] = sum-P[*b+*nb*i]; // compute and save p11(s,t) factor
	j = 0;
	surv = 1;
	getIndexI(S, index1, &UT[0], len, &j, &e); // determine first index
	while (j < e) { // loop through the sample until last index is reached
		n = *len-j; // count individuals at risk
		r = E[index1[j]]; // initialize event count
		d = 1-E[index1[j]]; // initialize censure count
		for (j++; j < e && S[index1[j]] == S[index1[j-1]]; j++) { // loop until time changes or last index is reached
			r += E[index1[j]]; // count the events
			d += 1-E[index1[j]]; // count the censures
		}
		if (n-r != 0) surv *= 1-(double)d/(n-r); // compute survival probability
	}
	getIndexI(S, index1, &UT[*nt-1], len, &j, &e); // determine index
	i = 0;
	while (j < e) { // loop through the sample until last index is reached
		n = *len-j; // count individuals at risk
		r = E[index1[j]]; // initialize event count
		d = 1-E[index1[j]]; // initialize censure count
		for (j++; j < e && S[index1[j]] == S[index1[j-1]]; j++) { // loop until time changes or last index is reached
			r += E[index1[j]]; // count the events
			d += 1-E[index1[j]]; // count the censures
		}
		if (n-r != 0) surv *= 1-(double)d/(n-r); // compute survival probability
		if (surv > 0) {
			k = *len-n;
			while (S[index1[k]] > UT[y]) y++;
			for (; k < j; k++) {
				if (T1[index1[k]] <= UT[0] && E[index1[k]]) for (i = 0; i < y; i++) P[*b+*nb*(i+*nt*3)] += 1/surv; // compute sum
				else if (E[index1[k]]) for (i = 0; i < y; i++) P[*b+*nb*(i+*nt)] += (T1[index1[k]] <= UT[i])/surv; // compute sum
			}
		}
	}
	while (j < *len) { // loop through the sample until last index is reached
		n = *len-j; // count individuals at risk
		r = E[index1[j]]; // initialize event count
		d = 1-E[index1[j]]; // initialize censure count
		for (j++; j < *len && S[index1[j]] == S[index1[j-1]]; j++) { // loop until time changes or last index is reached
			r += E[index1[j]]; // count the events
			d += 1-E[index1[j]]; // count the censures
		}
		if (n-r != 0) surv *= 1-(double)d/(n-r); // compute survival probability
		if (surv > 0) {
			k = *len-n;
			for (; k < j; k++) {
				if (T1[index1[k]] <= UT[0] && E[index1[k]]) for (i = 0; i < *nt; i++) P[*b+*nb*(i+*nt*3)] += 1/surv; // compute sum
				else if (E[index1[k]]) for (i = 0; i < *nt; i++) P[*b+*nb*(i+*nt)] += (T1[index1[k]] <= UT[i])/surv; // compute sum
			}
		}
	}
	for (i = *nt-1; i >= 0; i--) {
		P[*b+*nb*(i+*nt)] /= P[*b]; // compute and save p12(s,t)
		P[*b+*nb*i] /= P[*b]; // compute and save p11(s,t)
		P[*b+*nb*(i+*nt*2)] = 1-P[*b+*nb*i]-P[*b+*nb*(i+*nt)]; // compute and save p13(s,t)
		if (P[*b+*nb*(i+*nt*2)] < 0) {
			P[*b+*nb*(i+*nt)] = 1-P[*b+*nb*i]; // compute and save p12(s,t)
			P[*b+*nb*(i+*nt*2)] = 0; // save p13(s,t)
		}
		P[*b+*nb*(i+*nt*3)] /= P[*b+*nb*(*nt*3)]; // compute and save p22(s,t)
	}
} // transIPCW2I

/*
Author:
	Artur Agostinho Araujo <artur.stat@gmail.com>

Description:
	Computes a transition probability vector based
		on the inverse probability of censoring estimator.

Parameters:
	object			an object of class 'IPCW1'.
	UT				unique times vector.
	nboot			number of bootstrap samples.
	methodest		an integer indicating the desired method.

Return value:
	Returns a list where the first element is a
		(nboot)x(nt)x4 array of transition probabilities,
		and the second element is NULL.
*/

SEXP TransPROBIPCW1(
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
	PROTECT( P = alloc3DArray(REALSXP, *INTEGER(nboot), nt, 4) );
	PROTECT( list = NEW_LIST(2) );
	void (*func)(CintCP, CdoubleCP, CintCP, CdoubleCP, CintCP, CintCP, CintCP, CintCP, CdoubleCP, CintCP, doubleCP, CintCP);
	switch ( *INTEGER(methodest) ) {
		case 2:
			func = transIPCW2I;
			break;
		default:
			func = transIPCW1I;
	}
	#ifdef _OPENMP
	#pragma omp parallel if(*INTEGER(nboot) > 1) num_threads(global_num_threads)
	#endif
	{
		int b;
		int *index0 = (int*)malloc( len*sizeof(int) ); // allocate memory block
		int *index1 = (int*)malloc( len*sizeof(int) ); // allocate memory block
		double *WORK = (double*)malloc( len*sizeof(double) ); // allocate memory block
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
			order_d(REAL(T1), index0, len, FALSE, FALSE, WORK); // get permuation
			order_d(REAL(S), index1, len, FALSE, FALSE, WORK); // get permuation
			func(&len, REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), index0, index1, &nt, REAL(UT), INTEGER(nboot), REAL(P), &b); // compute transition probabilities
		}
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for (b = 1; b < *INTEGER(nboot); b++) {
			boot_ii(index0, index1, &len); // bootstrap indexes
			order_d(REAL(T1), index0, len, FALSE, FALSE, WORK); // get permuation
			order_d(REAL(S), index1, len, FALSE, FALSE, WORK); // get permuation
			func(&len, REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), index0, index1, &nt, REAL(UT), INTEGER(nboot), REAL(P), &b); // compute transition probabilities
		}
		free(index0); // free memory block
		free(index1); // free memory block
		free(WORK); // free memory block
	}
	SET_ELEMENT(list, 0, P);
	SET_ELEMENT(list, 1, R_NilValue);
	UNPROTECT(2);
	return list;
} // TransPROBIPCW1
