/*************************************/
/*** KMPW1 TRANSITION PROBABILITIES ***/
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
#include "logistic.h"

/*
Author:
	Artur Araujo <artur.stat@gmail.com>

Description:
	Computes the transition probabilities:
		p11(s,t) = P(Z>t|Z>s) = {1-P(Z<=t)}/{1-P(Z<=s)}
		p12(s,t) = P(Z<=t,T>t|Z>s) = {P(Z<=t)-P(Z<=s)-P(s<Z<=t,T<=t)}/{1-P(Z<=s)}
		p13(s,t) = 1-p11(s,t)-p12(s,t)
		p22(s,t) = P(Z<=t,T>t|Z<=s,T>s) = {P(Z<=s)-P(Z<=s,T<=t)}/{P(Z<=s)-P(T<=s)}

Parameters:
	len[in]			pointer to length of T1, M0, S and MS.
	T1[in]			pointer to T1 first element.
	M0[in]			pointer to M0 first element.
	S[in]			pointer to S first element.
	MS[in]			pointer to MS first element.
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
	Vectors T1, M0, S and MS must have the same length.
*/

static void transKMPW1I(
	CintCP len,
	Cdouble T1[*len],
	Cdouble M0[*len],
	Cdouble S[*len],
	Cdouble MS[*len],
	Cint index0[*len],
	Cint index1[*len],
	CintCP nt,
	Cdouble UT[*nt],
	CintCP nb,
	double P[*nb*(*nt)*4],
	CintCP b)
{
	register int i;
	int j = 0, e;
	double aux[7], p[2];
	getIndexI(T1, index0, &UT[0], len, &j, &e); // determine first index
	for (aux[0] = 1, aux[3] = 1; j < e; j++) { // loop through the sample until last index is reached
		aux[2] = M0[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[3] -= aux[2]; // compute probability
	}
	getIndexI(T1, index0, &UT[*nt-1], len, &j, &e); // determine last index
	for (i = 0, aux[4] = aux[3]; j < e; j++) { // loop through the sample until last index is reached
		p[0] = aux[4]/aux[3]; // compute transition probability
		p[1] = 1-p[0]; // compute probability
		while (T1[index0[j]] > UT[i]) {
			if (p[0] < 0) P[*b+*nb*i] = 0;
			else P[*b+*nb*i] = p[0]; // save p11(s,t)
			P[*b+*nb*(i+*nt)] = p[1]; // save probability
			P[*b+*nb*(i+*nt*3)] = 1;
			i++;
		}
		aux[2] = M0[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[4] -= aux[2]; // compute probability
	}
	p[0] = aux[4]/aux[3]; // compute and save transition probability
	p[1] = 1-p[0]; // compute and save probability
	if (p[0] < 0) p[0] = 0;
	for (; i < *nt; i++) {
		P[*b+*nb*i] = p[0]; // save p11(s,t)
		P[*b+*nb*(i+*nt)] = p[1]; // save probability
		P[*b+*nb*(i+*nt*3)] = 1;
	}
	j = 0;
	getIndexI(S, index1, &UT[0], len, &j, &e); // determine first index
	for (aux[0] = 1, aux[4] = 0; j < e; j++) { // loop through the sample until last index is reached
		aux[2] = MS[index1[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[4] += aux[2]; // compute probability
	}
	getIndexI(S, index1, &UT[*nt-1], len, &j, &e); // determine last index
	for (i = 0, aux[5] = 0, aux[6] = 0; j < e; j++) { // loop through the sample until last index is reached
		p[0] = aux[5]/aux[3]; // compute probability
		p[1] = aux[6]/(1-aux[3]-aux[4]); // compute probability
		while (S[index1[j]] > UT[i]) {
			P[*b+*nb*(i+*nt)] -= p[0]; // compute and save p12(s,t)
			if (P[*b+*nb*(i+*nt)] < 0) P[*b+*nb*(i+*nt)] = 0;
			P[*b+*nb*(i+*nt*2)] = 1-P[*b+*nb*i]-P[*b+*nb*(i+*nt)]; // compute and save p13(s,t)
			if (P[*b+*nb*(i+*nt*2)] < 0) {
				P[*b+*nb*(i+*nt)] = 1-P[*b+*nb*i]; // compute and save p12(s,t)
				P[*b+*nb*(i+*nt*2)] = 0; // compute and save p13(s,t)
			}
			P[*b+*nb*(i+*nt*3)] -= p[1]; // compute and save p22(s,t)
			if (P[*b+*nb*(i+*nt*3)] < 0) P[*b+*nb*(i+*nt*3)] = 0;
			i++;
		}
		aux[2] = MS[index1[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[5] += aux[2]*(T1[index1[j]] > UT[0]); // compute probability
		aux[6] += aux[2]*(T1[index1[j]] <= UT[0]); // compute probability
	}
	p[0] = aux[5]/aux[3]; // compute probability
	p[1] = aux[6]/(1-aux[3]-aux[4]); // compute probability
	for (; i < *nt; i++) {
		P[*b+*nb*(i+*nt)] -= p[0]; // compute and save p12(s,t)
		if (P[*b+*nb*(i+*nt)] < 0) P[*b+*nb*(i+*nt)] = 0;
		P[*b+*nb*(i+*nt*2)] = 1-P[*b+*nb*i]-P[*b+*nb*(i+*nt)]; // compute and save p13(s,t)
		if (P[*b+*nb*(i+*nt*2)] < 0) {
			P[*b+*nb*(i+*nt)] = 1-P[*b+*nb*i]; // compute and save p12(s,t)
			P[*b+*nb*(i+*nt*2)] = 0; // compute and save p13(s,t)
		}
		P[*b+*nb*(i+*nt*3)] -= p[1]; // compute and save p22(s,t)
		if (P[*b+*nb*(i+*nt*3)] < 0) P[*b+*nb*(i+*nt*3)] = 0;
	}
	return;
} // transKMPW1I

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
	len[in]			pointer to length of T1, M0, S and MS.
	T1[in]			pointer to T1 first element.
	M0[in]			pointer to M0 first element.
	S[in]			pointer to S first element.
	MS[in]			pointer to MS first element.
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
	Vectors T1, M0, S and MS must have the same length.
*/

static void transKMPW2I(
	CintCP len,
	Cdouble T1[*len],
	Cdouble M0[*len],
	Cdouble S[*len],
	Cdouble MS[*len],
	Cint index0[*len],
	Cint index1[*len],
	CintCP nt,
	Cdouble UT[*nt],
	CintCP nb,
	double P[*nb*(*nt)*4],
	CintCP b)
{
	register int i, k;
	int j = 0, e;
	double aux[4];
	getIndexI(T1, index0, &UT[0], len, &j, &e); // determine first index
	for (j = 0, aux[0] = 1, aux[3] = 0; j < e; j++) { // loop through the sample until last index is reached
		aux[2] = M0[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[3] += aux[2]; // compute probability
	}
	getIndexI(T1, index0, &UT[*nt-1], len, &j, &e); // determine last index
	for (i = 0; j < e; j++) { // loop through the sample until last index is reached
		while (T1[index0[j]] > UT[i]) {
			P[*b+*nb*i] = aux[3];
			P[*b+*nb*(i+*nt)] = 0;
			P[*b+*nb*(i+*nt*3)] = 0;
			i++;
		}
		aux[2] = M0[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[3] += aux[2]; // compute probability
	}
	for (; i < *nt; i++) {
		P[*b+*nb*i] = aux[3];
		P[*b+*nb*(i+*nt)] = 0;
		P[*b+*nb*(i+*nt*3)] = 0;
	}
	for (; j < *len; j++) { // loop through the sample until last index is reached
		aux[2] = M0[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[3] += aux[2]; // compute probability
	}
	for (i = 0; i < *nt; i++) P[*b+*nb*i] = aux[3]-P[*b+*nb*i];
	j = 0;
	getIndexI(S, index1, &UT[0], len, &j, &e); // determine first index
	for (aux[0] = 1; j < e; j++) { // loop through the sample until last index is reached
		aux[1] = 1-MS[index1[j]]/(*len-j); // compute needed factor
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
	}
	getIndexI(S, index1, &UT[*nt-1], len, &j, &e); // determine first index
	for (k = 0; j < e; j++) { // loop through the sample until last index is reached
		aux[2] = MS[index1[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		while (S[index1[j]] > UT[k]) k++;
		if (T1[index1[j]] <= UT[0]) for (i = 0; i < k; i++) P[*b+*nb*(i+*nt*3)] += aux[2];
		else for (i = 0; i < k; i++) P[*b+*nb*(i+*nt)] += aux[2]*(T1[index1[j]] <= UT[i]);
	}
	for (; j < *len; j++) { // loop through the sample until last index is reached
		aux[2] = MS[index1[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		if (T1[index1[j]] <= UT[0]) for (i = 0; i < *nt; i++) P[*b+*nb*(i+*nt*3)] += aux[2];
		else for (i = 0; i < *nt; i++) P[*b+*nb*(i+*nt)] += aux[2]*(T1[index1[j]] <= UT[i]);
	}
	for (i = *nt-1; i >= 0; i--) {
		P[*b+*nb*(i+*nt)] /= P[*b]; // compute and save p12(s,t)
		P[*b+*nb*i] /= P[*b]; // compute and save p11(s,t)
		P[*b+*nb*(i+*nt*2)] = 1-P[*b+*nb*i]-P[*b+*nb*(i+*nt)]; // compute and save p13(s,t)
		if (P[*b+*nb*(i+*nt*2)] < 0) {
			P[*b+*nb*(i+*nt)] = 1-P[*b+*nb*i]; // compute and save p12(s,t)
			P[*b+*nb*(i+*nt*2)] = 0; // compute and save p13(s,t)
		}
		P[*b+*nb*(i+*nt*3)] /= P[*b+*nb*(*nt*3)]; // compute and save p22(s,t)
	}
	return;
} // transKMPW2I

/*
Author:
	Artur Araujo <artur.stat@gmail.com>

Description:
	Computes a transition probability vector based
		on the presmoothed Kaplan-Meier weights estimator.

Parameters:
	object			an object of class 'KMPW1'.
	UT			unique times vector.
	nboot			number of bootstrap samples.
	methodest		an integer indicating the desired method.

Return value:
	Returns a list where the first element is a
		(nboot)x(nt)x4 array of transition probabilities,
		and the second element is NULL.
*/

SEXP TransPROBKMPW1(
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
	int len = GET_LENGTH(T1), nt = GET_LENGTH(UT), t, nth = 1;
	SEXP P, list;
	PROTECT( P = alloc3DArray(REALSXP, *INTEGER(nboot), nt, 4) );
	PROTECT( list = NEW_LIST(2) );
	void (*func)(CintCP, CdoubleCP, CdoubleCP, CdoubleCP, CdoubleCP, CintCP, CintCP, CintCP, CdoubleCP, CintCP, doubleCP, CintCP);
	switch ( *INTEGER(methodest) ) {
		case 2:
			func = transKMPW2I;
			break;
		default:
			func = transKMPW1I;
	}
	int b, n0 = 2, ns = 3, maxit = 30;
	double *J = (double*)malloc( len*sizeof(double) ); // allocate memory block
	if (J == NULL) error("TransPROBKMPW1: No more memory\n");
	double *X0[2] = {J, REAL(T1)};
	double *XS[3] = {J, REAL(T1), REAL(S)};
	double epsilon = 1e-8;
	for (b = 0; b < len; b++) J[b] = 1; // initialize J vector
	if (*INTEGER(nboot) > 1) nth = global_num_threads;
	int **index0 = (int**)malloc( nth*sizeof(int*) ); // allocate memory block
	if (index0 == NULL) error("TransPROBKMPW1: No more memory\n");
	int **index1 = (int**)malloc( nth*sizeof(int*) ); // allocate memory block
	if (index1 == NULL) error("TransPROBKMPW1: No more memory\n");
	double **M0 = (double**)malloc( nth*sizeof(double*) ); // allocate memory block
	if (M0 == NULL) error("TransPROBKMPW1: No more memory\n");
	double **MS = (double**)malloc( nth*sizeof(double*) ); // allocate memory block
	if (MS == NULL) error("TransPROBKMPW1: No more memory\n");
	double **WORK0 = (double**)malloc( nth*sizeof(double*) ); // allocate memory block
	if (WORK0 == NULL) error("TransPROBKMPW1: No more memory\n");
	double **WORK1 = (double**)malloc( nth*sizeof(double*) ); // allocate memory block
	if (WORK1 == NULL) error("TransPROBKMPW1: No more memory\n");
	logitW **WORK = (logitW**)malloc( nth*sizeof(logitW*) ); // allocate memory block
	if (WORK == NULL) error("TransPROBKMPW1: No more memory\n");
	for (t = 0; t < nth; t++) { // allocate per thread memory
		if ( ( index0[t] = (int*)malloc( len*sizeof(int) ) ) == NULL ) error("TransPROBKMPW1: No more memory\n");
		if ( ( index1[t] = (int*)malloc( len*sizeof(int) ) ) == NULL ) error("TransPROBKMPW1: No more memory\n");
		if ( ( M0[t] = (double*)malloc( len*sizeof(double) ) ) == NULL ) error("TransPROBKMPW1: No more memory\n");
		if ( ( MS[t] = (double*)malloc( len*sizeof(double) ) ) == NULL ) error("TransPROBKMPW1: No more memory\n");
		if ( ( WORK0[t] = (double*)malloc( len*sizeof(double) ) ) == NULL ) error("TransPROBKMPW1: No more memory\n");
		if ( ( WORK1[t] = (double*)malloc( len*sizeof(double) ) ) == NULL ) error("TransPROBKMPW1: No more memory\n");
		WORK[t] = logitW_Create(&ns);
	}
	#ifdef _OPENMP
	#pragma omp parallel num_threads(nth) private(b, t)
	#endif
	{
		int conv;
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
			predict_logit(&len, index0[t], INTEGER(E1), M0[t], &n0, X0, &maxit, &epsilon, &conv, WORK[t]); // compute M0 predicted values
			predict_logit(&len, index1[t], INTEGER(E), MS[t], &ns, XS, &maxit, &epsilon, &conv, WORK[t]); // compute MS predicted values
			order_dd(REAL(T1), M0[t], index0[t], len, FALSE, FALSE, TRUE, WORK0[t], WORK1[t]); // get permuation
			order_dd(REAL(S), MS[t], index1[t], len, FALSE, FALSE, TRUE, WORK0[t], WORK1[t]); // get permuation
			func(&len, REAL(T1), M0[t], REAL(S), MS[t], index0[t], index1[t], &nt, REAL(UT), INTEGER(nboot), REAL(P), &b); // compute transition probabilities
		}
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for (b = 1; b < *INTEGER(nboot); b++) {
			boot_ii(RngArray[t], &len, index0[t], index1[t]); // bootstrap indexes
			predict_logit(&len, index0[t], INTEGER(E1), M0[t], &n0, X0, &maxit, &epsilon, &conv, WORK[t]); // compute M0 predicted values
			predict_logit(&len, index1[t], INTEGER(E), MS[t], &ns, XS, &maxit, &epsilon, &conv, WORK[t]); // compute MS predicted values
			order_dd(REAL(T1), M0[t], index0[t], len, FALSE, FALSE, TRUE, WORK0[t], WORK1[t]); // get permuation
			order_dd(REAL(S), MS[t], index1[t], len, FALSE, FALSE, TRUE, WORK0[t], WORK1[t]); // get permuation
			func(&len, REAL(T1), M0[t], REAL(S), MS[t], index0[t], index1[t], &nt, REAL(UT), INTEGER(nboot), REAL(P), &b); // compute transition probabilities
		}
	}
	for (t = nth-1; t >= 0; t--) {
		free(index0[t]); // free memory block
		free(index1[t]); // free memory block
		free(M0[t]); // free memory block
		free(MS[t]); // free memory block
		free(WORK0[t]); // free memory block
		free(WORK1[t]); // free memory block
		logitW_Delete(WORK[t]);
	}
	free(index0); // free memory block
	free(index1); // free memory block
	free(M0); // free memory block
	free(MS); // free memory block
	free(WORK0); // free memory block
	free(WORK1); // free memory block
	free(WORK); // free memory block
	free(J); // free memory block
	SET_ELEMENT(list, 0, P);
	SET_ELEMENT(list, 1, R_NilValue);
	UNPROTECT(2);
	return list;
} // TransPROBKMPW1
