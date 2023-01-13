/*************************************/
/*** KMPW2 TRANSITION PROBABILITIES ***/
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
    p12(s,t) = 1-p11(s,t)-p13(s,t)
    p13(s,t) = P(T<=t|Z>s) = P(Z>s,T<=t)/{1-P(Z<=s)}
    p23(s,t) = P(T<=t|Z<=s,T>s) = P(Z<=s,s<T<=t)/{P(Z<=s)-P(T<=s)}

Parameters:
  len[in]           pointer to length of T1, M0, S and M.
  T1[in]            pointer to T1 first element.
  M0[in]            pointer to M0 first element.
  S[in]             pointer to S first element.
  M[in]             pointer to M first element.
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
  Vectors T1, M0, S and M must have the same length.
*/

static void transKMPW3I(
	CintCP len,
	Cdouble T1[*len],
	Cdouble M0[*len],
	Cdouble S[*len],
	Cdouble M[*len],
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
	double aux[4], p[2];
	getIndexI(T1, index0, &UT[0], len, &j, &e); // determine first index
	for (aux[0] = 1, aux[3] = 1; j < e; j++) { // loop through the sample until last index is reached
		aux[2] = M0[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[3] -= aux[2]; // compute probability
	}
	getIndexI(T1, index0, &UT[*nt-1], len, &j, &e); // determine last index
	for (i = 0; j < e; j++) { // loop through the sample until last index is reached
		while (T1[index0[j]] > UT[i]) P[*b+*nb*i++] = aux[3];
		aux[2] = M0[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[3] -= aux[2]; // compute probability
	}
	for (; i < *nt; i++) P[*b+*nb*i] = aux[3];
	j = 0;
	getIndexI(S, index1, &UT[0], len, &j, &e); // determine first index
	for (aux[0] = 1, aux[3] = 0; j < e; j++) { // loop through the sample until last index is reached
		aux[2] = M[index1[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[3] += aux[2]; // compute probability
	}
	getIndexI(S, index1, &UT[*nt-1], len, &j, &e); // determine last index
	for (i = 0, p[0] = 0, p[1] = 0; j < e; j++) { // loop through the sample until last index is reached
		while (S[index1[j]] > UT[i]) {
			P[*b+*nb*(i+*nt*2)] = p[0];
			P[*b+*nb*(i+*nt*3)] = p[1];
			i++;
		}
		aux[2] = M[index1[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		if (T1[index1[j]] <= UT[0]) p[1] += aux[2]; // compute probability
		else p[0] += aux[2]; // compute probability
	}
	for (; i < *nt; i++) {
		P[*b+*nb*(i+*nt*2)] = p[0];
		P[*b+*nb*(i+*nt*3)] = p[1];
	}
	for (i = *nt-1; i >= 0; i--) {
		P[*b+*nb*(i+*nt*2)] /= P[*b]; // compute and save p13(s,t)
		P[*b+*nb*(i+*nt*3)] /= 1-P[*b]-aux[3]; // compute and save p23(s,t)
		if (P[*b+*nb*(i+*nt*3)] > 1) P[*b+*nb*(i+*nt*3)] = 1;
		P[*b+*nb*i] /= P[*b]; // compute and save p11(s,t)
		P[*b+*nb*(i+*nt)] = 1-P[*b+*nb*i]-P[*b+*nb*(i+*nt*2)]; // compute and save p12(s,t)
		if (P[*b+*nb*(i+*nt)] < 0) {
			P[*b+*nb*(i+*nt*2)] = 1-P[*b+*nb*i]; // compute and save p13(s,t)
			P[*b+*nb*(i+*nt)] = 0; // compute and save p12(s,t)
		}
	}
	return;
} // transKMPW3I

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Computes the transition probabilities:
    p11(s,t) = P(Z>t|Z>s) = P(Z>t)/P(Z>s)
    p12(s,t) = 1-p11(s,t)-p13(s,t)
    p13(s,t) = P(T<=t|Z>s) = P(Z>s,T<=t)/P(Z>s)
    p23(s,t) = P(T<=t|Z<=s,T>s) = P(Z<=s,s<T<=t)/P(Z<=s,T>s)

Parameters:
  len[in]           pointer to length of T1, M0, S and M.
  T1[in]            pointer to T1 first element.
  M0[in]            pointer to M0 first element.
  S[in]             pointer to S first element.
  M[in]             pointer to M first element.
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
  Vectors T1, M0, S and M must have the same length.
*/

static void transKMPW4I(
	CintCP len,
	Cdouble T1[*len],
	Cdouble M0[*len],
	Cdouble S[*len],
	Cdouble M[*len],
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
	double aux[4], p[2];
	getIndexI(T1, index0, &UT[0], len, &j, &e); // determine first index
	for (aux[0] = 1, aux[3] = 0; j < e; j++) { // loop through the sample until last index is reached
		aux[2] = M0[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[3] += aux[2]; // compute probability
	}
	getIndexI(T1, index0, &UT[*nt-1], len, &j, &e); // determine last index
	for (i = 0; j < e; j++) { // loop through the sample until last index is reached
		while (T1[index0[j]] > UT[i]) P[*b+*nb*i++] = aux[3];
		aux[2] = M0[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		aux[3] += aux[2]; // compute probability
	}
	for (; i < *nt; i++) P[*b+*nb*i] = aux[3];
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
		aux[1] = 1-M[index1[j]]/(*len-j); // compute needed factor
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
	}
	getIndexI(S, index1, &UT[*nt-1], len, &j, &e); // determine last index
	for (i = 0, p[0] = 0, p[1] = 0; j < e; j++) { // loop through the sample until last index is reached
		while (S[index1[j]] > UT[i]) {
			P[*b+*nb*(i+*nt*2)] = p[0]; // save probability
			P[*b+*nb*(i+*nt*3)] = p[1]; // save probability
			i++;
		}
		aux[2] = M[index1[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		if (T1[index1[j]] <= UT[0]) p[1] += aux[2]; // compute probability
		else p[0] += aux[2]; // compute probability
	}
	for (; i < *nt; i++) {
		P[*b+*nb*(i+*nt*2)] = p[0];
		P[*b+*nb*(i+*nt*3)] = p[1];
	}
	for (; j < *len; j++) { // loop through the sample until last index is reached
		aux[2] = M[index1[j]]/(*len-j); // compute needed factor
		aux[1] = 1-aux[2]; // factor needed for the computation
		aux[2] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		if (T1[index1[j]] <= UT[0]) p[1] += aux[2]; // compute probability
	}
	for (i = *nt-1; i >= 0; i--) {
		P[*b+*nb*(i+*nt*2)] /= P[*b]; // compute and save p13(s,t)
		P[*b+*nb*i] /= P[*b]; // compute and save p11(s,t)
		P[*b+*nb*(i+*nt)] = 1-P[*b+*nb*i]-P[*b+*nb*(i+*nt*2)]; // compute and save p12(s,t)
		if (P[*b+*nb*(i+*nt)] < 0) {
			P[*b+*nb*(i+*nt*2)] = 1-P[*b+*nb*i]; // compute and save p13(s,t)
			P[*b+*nb*(i+*nt)] = 0; // compute and save p12(s,t)
		}
		P[*b+*nb*(i+*nt*3)] /= p[1]; // compute and save p23(s,t)
	}
	return;
} // transKMPW4I

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Computes a transition probability vector based
    on the presmoothed Kaplan-Meier weights estimator.

Parameters:
  object            an object of class 'KMPW2'.
  UT                unique times vector.
  nboot             number of bootstrap samples.
  methodest         an integer indicating the desired method.

Return value:
  Returns a list where the first element is a
    (nboot)x(nt)x4 array of transition probabilities,
    and the second element is NULL.
*/

SEXP TransPROBKMPW2(
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
		case 3:
			func = transKMPW3I;
			break;
		default:
			func = transKMPW4I;
	}
	register int j;
	int n0 = 2, n1 = 3, maxit = 30;
	double *J = (double*)malloc( (unsigned int)len*sizeof(double) ); // allocate memory block
	if (J == NULL) error("TransPROBKMPW2: No more memory\n");
	double *X0[2] = {J, REAL(T1)};
	double *XM[3] = {J, REAL(T1), REAL(S)};
	double epsilon = 1e-8;
	for (j = 0; j < len; j++) J[j] = 1; // initialize J vector
	if (*INTEGER(nboot) > 1) nth = global_num_threads;
	int **index0 = (int**)malloc( (unsigned int)nth*sizeof(int*) ); // allocate memory block
	if (index0 == NULL) error("TransPROBKMPW2: No more memory\n");
	int **index1 = (int**)malloc( (unsigned int)nth*sizeof(int*) ); // allocate memory block
	if (index1 == NULL) error("TransPROBKMPW2: No more memory\n");
	double **M0 = (double**)malloc( (unsigned int)nth*sizeof(double*) ); // allocate memory block
	if (M0 == NULL) error("TransPROBKMPW2: No more memory\n");
	double **M = (double**)malloc( (unsigned int)nth*sizeof(double*) ); // allocate memory block
	if (M == NULL) error("TransPROBKMPW2: No more memory\n");
	int **subset0 = (int**)malloc( (unsigned int)nth*sizeof(int*) ); // allocate memory block
	if (subset0 == NULL) error("TransPROBKMPW2: No more memory\n");
	int **subset1 = (int**)malloc( (unsigned int)nth*sizeof(int*) ); // allocate memory block
	if (subset1 == NULL) error("TransPROBKMPW2: No more memory\n");
	double **WORK0 = (double**)malloc( (unsigned int)nth*sizeof(double*) ); // allocate memory block
	if (WORK0 == NULL) error("TransPROBKMPW2: No more memory\n");
	double **WORK1 = (double**)malloc( (unsigned int)nth*sizeof(double*) ); // allocate memory block
	if (WORK1 == NULL) error("TransPROBKMPW2: No more memory\n");
	logitW **WORK = (logitW**)malloc( (unsigned int)nth*sizeof(logitW*) ); // allocate memory block
	if (WORK == NULL) error("TransPROBKMPW2: No more memory\n");
	for (t = 0; t < nth; t++) { // allocate per thread memory
		if ( ( index0[t] = (int*)malloc( (unsigned int)len*sizeof(int) ) ) == NULL ) error("TransPROBKMPW2: No more memory\n");
		if ( ( index1[t] = (int*)malloc( (unsigned int)len*sizeof(int) ) ) == NULL ) error("TransPROBKMPW2: No more memory\n");
		if ( ( M0[t] = (double*)malloc( (unsigned int)len*sizeof(double) ) ) == NULL ) error("TransPROBKMPW2: No more memory\n");
		if ( ( M[t] = (double*)malloc( (unsigned int)len*sizeof(double) ) ) == NULL ) error("TransPROBKMPW2: No more memory\n");
		if ( ( subset0[t] = (int*)malloc( (unsigned int)len*sizeof(int) ) ) == NULL ) error("TransPROBKMPW2: No more memory\n");
		if ( ( subset1[t] = (int*)malloc( (unsigned int)len*sizeof(int) ) ) == NULL ) error("TransPROBKMPW2: No more memory\n");
		if ( ( WORK0[t] = (double*)malloc( (unsigned int)len*sizeof(double) ) ) == NULL ) error("TransPROBKMPW2: No more memory\n");
		if ( ( WORK1[t] = (double*)malloc( (unsigned int)len*sizeof(double) ) ) == NULL ) error("TransPROBKMPW2: No more memory\n");
		WORK[t] = logitW_Create(&n1);
	}
	#ifdef _OPENMP
	#pragma omp parallel num_threads(nth) private(j, t)
	#endif
	{
		int b, s0, s1, conv;
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
			for (j = 0, s0 = 0, s1 = 0; j < len; j++) { // subset sample
				if (REAL(T1)[index1[t][j]] == REAL(S)[index1[t][j]]) subset0[t][s0++] = index1[t][j];
				else subset1[t][s1++] = index1[t][j];
			}
			predict_logit(&s0, subset0[t], INTEGER(E1), M[t], &n0, X0, &maxit, &epsilon, &conv, WORK[t]); // compute M predicted values
			predict_logit(&s1, subset1[t], INTEGER(E), M[t], &n1, XM, &maxit, &epsilon, &conv, WORK[t]); // compute M predicted values
			order_dd(REAL(T1), M0[t], index0[t], len, FALSE, FALSE, TRUE, WORK0[t], WORK1[t]); // get permuation
			order_dd(REAL(S), M[t], index1[t], len, FALSE, FALSE, TRUE, WORK0[t], WORK1[t]); // get permuation
			func(&len, REAL(T1), M0[t], REAL(S), M[t], index0[t], index1[t], &nt, REAL(UT), INTEGER(nboot), REAL(P), &b); // compute transition probabilities
		}
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for (b = 1; b < *INTEGER(nboot); b++) {
			boot_ii(RngArray[t], &len, index0[t], index1[t]); // bootstrap indexes
			predict_logit(&len, index0[t], INTEGER(E1), M0[t], &n0, X0, &maxit, &epsilon, &conv, WORK[t]); // compute M0 predicted values
			for (j = 0, s0 = 0, s1 = 0; j < len; j++) { // subset sample
				if (REAL(T1)[index1[t][j]] == REAL(S)[index1[t][j]]) subset0[t][s0++] = index1[t][j];
				else subset1[t][s1++] = index1[t][j];
			}
			predict_logit(&s0, subset0[t], INTEGER(E1), M[t], &n0, X0, &maxit, &epsilon, &conv, WORK[t]); // compute M predicted values
			predict_logit(&s1, subset1[t], INTEGER(E), M[t], &n1, XM, &maxit, &epsilon, &conv, WORK[t]); // compute M predicted values
			order_dd(REAL(T1), M0[t], index0[t], len, FALSE, FALSE, TRUE, WORK0[t], WORK1[t]); // get permuation
			order_dd(REAL(S), M[t], index1[t], len, FALSE, FALSE, TRUE, WORK0[t], WORK1[t]); // get permuation
			func(&len, REAL(T1), M0[t], REAL(S), M[t], index0[t], index1[t], &nt, REAL(UT), INTEGER(nboot), REAL(P), &b); // compute transition probabilities
		}
	}
	for (t = nth-1; t >= 0; t--) {
		free(index0[t]); // free memory block
		free(index1[t]); // free memory block
		free(M0[t]); // free memory block
		free(M[t]); // free memory block
		free(subset0[t]); // free memory block
		free(subset1[t]); // free memory block
		free(WORK0[t]); // free memory block
		free(WORK1[t]); // free memory block
		logitW_Delete(WORK[t]);
	}
	free(index0); // free memory block
	free(index1); // free memory block
	free(M0); // free memory block
	free(M); // free memory block
	free(subset0); // free memory block
	free(subset1); // free memory block
	free(WORK0); // free memory block
	free(WORK1); // free memory block
	free(WORK); // free memory block
	free(J); // free memory block
	SET_ELEMENT(list, 0, P);
	SET_ELEMENT(list, 1, R_NilValue);
	UNPROTECT(2);
	return list;
} // TransPROBKMPW2
