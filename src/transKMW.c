/************************************/
/*** KMW TRANSITION PROBABILITIES ***/
/************************************/

#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdlib.h>
#include <time.h>
#include <Rdefines.h>
#include "defines.h"
#include "boot.h"
#include "rthreads.h"
#include "sort.h"
#include "transKMW1.h"
#include "transKMW2.h"

/*
Author:
	Artur Agostinho Araujo <artur.stat@gmail.com>

Description:
	Computes a transition probability vector based
		on the Kaplan-Meier weights estimator.

Parameters:
	object			an object of class 'KMW1' or class 'KMW2'.
	UT				unique times vector.
	nboot			number of bootstrap samples.
	methodest		an integer indicating the desired method.

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
	int len = GET_LENGTH(T1), nt = GET_LENGTH(UT);
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
	#ifdef _OPENMP
	#pragma omp parallel if(*INTEGER(nboot) > 1) num_threads(global_num_threads)
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
} // TransPROBKMW
