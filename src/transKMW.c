/************************************/
/*** KMW TRANSITION PROBABILITIES ***/
/************************************/

#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdlib.h>
#include <Rdefines.h>
#include "defines.h"
#include "RngStream.h"
#include "RngArray.h"
#include "RngBoot.h"
#include "rthreads.h"
#include "sort.h"
#include "transKMW1.h"
#include "transKMW2.h"

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Computes a transition probability vector based
    on the Kaplan-Meier weights estimator.

Parameters:
  object            an object of class 'KMW1' or class 'KMW2'.
  UT                unique times vector.
  nboot             number of bootstrap samples.
  methodest         an integer indicating the desired method.

Return value:
  Returns a list where the first element is a
    (nboot)x(nt)x4 array of transition probabilities,
    and the second element is NULL.
*/

SEXP TransPROBKMW(
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
	void (*func)(CintCP, CdoubleCP, CintCP, CdoubleCP, CintCP, CintCP, CintCP, CintCP, CdoubleCP, CintCP, doubleCP, CintCP);
	switch ( *INTEGER(methodest) ) {
		case 2:
			func = transKMW2I;
			break;
		case 3:
			func = transKMW3I;
			break;
		case 4:
			func = transKMW4I;
			break;
		default:
			func = transKMW1I;
	}
	if (*INTEGER(nboot) > 1) nth = global_num_threads;
	int **index0 = (int**)malloc( (unsigned int)nth*sizeof(int*) ); // allocate memory block
	if (index0 == NULL) error("TransPROBKMW: No more memory\n");
	int **index1 = (int**)malloc( (unsigned int)nth*sizeof(int*) ); // allocate memory block
	if (index1 == NULL) error("TransPROBKMW: No more memory\n");
	double **WORK0 = (double**)malloc( (unsigned int)nth*sizeof(double*) ); // allocate memory block
	if (WORK0 == NULL) error("TransPROBKMW: No more memory\n");
	int **WORK1 = (int**)malloc( (unsigned int)nth*sizeof(int*) ); // allocate memory block
	if (WORK1 == NULL) error("TransPROBKMW: No more memory\n");
	for (t = 0; t < nth; t++) { // allocate per thread memory
		if ( ( index0[t] = (int*)malloc( (unsigned int)len*sizeof(int) ) ) == NULL ) error("TransPROBKMW: No more memory\n");
		if ( ( index1[t] = (int*)malloc( (unsigned int)len*sizeof(int) ) ) == NULL ) error("TransPROBKMW: No more memory\n");
		if ( ( WORK0[t] = (double*)malloc( (unsigned int)len*sizeof(double) ) ) == NULL ) error("TransPROBKMW: No more memory\n");
		if ( ( WORK1[t] = (int*)malloc( (unsigned int)len*sizeof(int) ) ) == NULL ) error("TransPROBKMW: No more memory\n");
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
			order_di(REAL(T1), INTEGER(E1), index0[t], len, FALSE, FALSE, TRUE, WORK0[t], WORK1[t]); // get permuation
			order_di(REAL(S), INTEGER(E), index1[t], len, FALSE, FALSE, TRUE, WORK0[t], WORK1[t]); // get permuation
			func(&len, REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), index0[t], index1[t], &nt, REAL(UT), INTEGER(nboot), REAL(P), &b); // compute transition probabilities
		}
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for (b = 1; b < *INTEGER(nboot); b++) {
			boot_ii(RngArray[t], &len, index0[t], index1[t]); // bootstrap indexes
			order_di(REAL(T1), INTEGER(E1), index0[t], len, FALSE, FALSE, TRUE, WORK0[t], WORK1[t]); // get permuation
			order_di(REAL(S), INTEGER(E), index1[t], len, FALSE, FALSE, TRUE, WORK0[t], WORK1[t]); // get permuation
			func(&len, REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), index0[t], index1[t], &nt, REAL(UT), INTEGER(nboot), REAL(P), &b); // compute transition probabilities
		}
	}
	for (t = nth-1; t >= 0; t--) {
		free(index0[t]); // free memory block
		free(index1[t]); // free memory block
		free(WORK0[t]); // free memory block
		free(WORK1[t]); // free memory block
	}
	free(index0); // free memory block
	free(index1); // free memory block
	free(WORK0); // free memory block
	free(WORK1); // free memory block
	SET_ELEMENT(list, 0, P);
	SET_ELEMENT(list, 1, R_NilValue);
	UNPROTECT(2);
	return list;
} // TransPROBKMW
