/*************************************/
/*** LIN2 TRANSITION PROBABILITIES ***/
/*************************************/

#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <Rdefines.h>
#include "defines.h"
#include "boot.h"
#include "get.h"
#include "rthreads.h"
#include "sort.h"
#include "wikmsurv.h"
#include "wtypefunc.h"
#include "window.h"

/*
Author:
	Artur Agostinho Araujo <artur.stat@gmail.com>

Description:
	Computes the conditional transition probabilities:
		p11(s,t|X) = P(Z>t|Z>s|X) = P(Z>t|X)/P(Z>s|X)
		p12(s,t|X) = P(Z<=t,T>t|Z>s|X) = P(s<Z<=t,T>t|X)/P(Z>s|X)
		p13(s,t|X) = 1-p11(s,t|X)-p12(s,t|X)
		p22(s,t|X) = P(Z<=t,T>t|Z<=s,T>s|X) = P(Z<=s,T>t|X)/P(Z<=s,T>s|X)

Parameters:
	len[in]			pointer to length of T1, E1, S, E, X and SW->ptr.
	T1[in]			pointer to T1 first element.
	E1[in]			pointer to E1 first element.
	S[in]			pointer to S first element.
	E[in]			pointer to E first element.
	X[in]			pointer to X first element.
	SW[in]			pointer to a weights stype structure.
	index0[in]		pointer to index0 first element.
	index1[in]		pointer to index1 first element.
	nt[in]			pointer to length of UT and number of columns of P.
	UT[in]			pointer to unique times vector.
	nx[in]			pointer to length of UX and number of faces of P.
	UX[in]			pointer to unique covariate vector.
	h[in]			pointer to bandwidth parameter.
	kfunc[in]		pointer to kernel density function.
	wfunc[in]		pointer to weights function.
	nb[in]			pointer to number of rows of P.
	P[out]			pointer to a (nb)x(nt)x(nx)x4 probability array.
	b[in]			pointer to row index.

Return value:
	This function doesn't return a value.

Remarks:
	Vector index0 must indicate the permutation of vector T1
		sorted by ascending order.
	Vector index1 must indicate the permutation of vector S
		sorted by ascending order.
	Vectors T1, E1, S, E, X and SW->ptr must have the same length.
*/

static void transLIN2I(
	CintCP len,
	Cdouble T1[*len],
	Cint E1[*len],
	Cdouble S[*len],
	Cint E[*len],
	Cdouble X[*len],
	CstypeCP SW,
	Cint index0[*len],
	Cint index1[*len],
	CintCP nt,
	Cdouble UT[*nt],
	CintCP nx,
	Cdouble UX[*nx],
	CdoubleCP h,
	Kfunc kfunc,
	Wfunc wfunc,
	CintCP nb,
	double P[*nb*(*nt)*(*nx)*4],
	CintCP b)
{
	const int64_t nbt = *nb*(*nt), nbtx = nbt*(*nx);
	int z, e[4];
	z = 0;
	getIndexI(T1, index0, &UT[0], len, &z, &e[0]); // determine first index
	z = e[0];
	getIndexI(T1, index0, &UT[*nt-1], len, &z, &e[1]); // determine last index
	z = 0;
	getIndexI(S, index1, &UT[0], len, &z, &e[2]); // determine first index
	z = e[2];
	getIndexI(S, index1, &UT[*nt-1], len, &z, &e[3]); // determine last index
	#ifdef _OPENMP
	#pragma omp parallel if(*b < 1) num_threads(global_num_threads) private(z)
	#endif
	{
		register int i, j;
		register int64_t k;
		int64_t k0;
		double p, sum[3];
		double *K = (double*)malloc( *len*sizeof(double) ); // allocate memory block
		double *SV = (double*)malloc( *len*sizeof(double) ); // allocate memory block
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for (j = 0; j < *nx; j++) {
			wfunc(X, SW, index0, &UX[j], h, K, kfunc); // compute weights
			wikmsurv(len, T1, E1, K, index0, &e[1], SV); // compute conditional survival probabilities
			for (z = 0, sum[0] = 0; z < e[0]; z++) sum[0] += K[index0[z]];
			for (sum[2] = sum[0]; z < *len; z++) sum[2] += K[index0[z]];
			z = e[0]; i = 0; k0 = *b+nbt*j; k = k0; sum[1] = 0;
			if (z == 0) {
				while (T1[index0[0]] > UT[i]) {
					P[k] = sum[2]; // k = *b+*nb*(i+*nt*j) = *b+nbt*j+*nb*i
					P[k+nbtx] = 0; // k+nbtx = *b+*nb*( i+*nt*j+*nt*(*nx) ) = *b+nbt*j+*nb*i+nbtx
					i++;
					k += *nb;
				}
				sum[1] += K[index0[0]];
				z++;
			}
			while (z < e[1]) {
				p = (sum[2]-sum[0]-sum[1])/SV[index0[z-1]];
				while (T1[index0[z]] > UT[i]) {
					P[k] = p; // k = *b+*nb*(i+*nt*j) = *b+nbt*j+*nb*i
					P[k+nbtx] = sum[1]; // k+nbtx = *b+*nb*( i+*nt*j+*nt*(*nx) ) = *b+nbt*j+*nb*i+nbtx
					i++;
					k += *nb;
				}
				sum[1] += K[index0[z]];
				z++;
			}
			p = (sum[2]-sum[0]-sum[1])/SV[index0[z-1]];
			while (k < k0+nbt) { // needed for bootstrap
				P[k] = p; // k = *b+*nb*(i+*nt*j) = *b+nbt*j+*nb*i
				P[k+nbtx] = sum[1]; // k+nbtx = *b+*nb*( i+*nt*j+*nt*(*nx) ) = *b+nbt*j+*nb*i+nbtx
				k += *nb;
			}
			wikmsurv(len, S, E, K, index1, &e[3], SV); // compute conditional survival probabilities
			for (z = 0, sum[1] = 0; z < e[2]; z++) sum[1] += K[index1[z]];
			i = 0; k = k0; sum[0] -= sum[1]; sum[1] = 0;
			if (z == 0) {
				while (S[index1[0]] > UT[i]) {
					P[k+nbtx*3] = sum[0]; // k+nbtx*3 = *b+*nb*(i+*nt*j+*nt*(*nx)*3) = *b+nbt*j+*nb*i+nbtx*3
					i++;
					k += *nb;
				}
				if (T1[index1[0]] <= UT[0]) sum[0] -= K[index1[0]]; // compute sum
				else sum[1] += K[index1[0]]; // compute sum
				z++;
			}
			while (z < e[3]) {
				p = sum[0]/SV[index1[z-1]];
				while (S[index1[z]] > UT[i]) {
					P[k+nbtx] -= sum[1]; // k+nbtx = *b+*nb*( i+*nt*j+*nt*(*nx) ) = *b+nbt*j+*nb*i+nbtx
					P[k+nbtx] /= SV[index1[z-1]];
					P[k+nbtx*3] = p; // k+nbtx*3 = *b+*nb*(i+*nt*j+*nt*(*nx)*3) = *b+nbt*j+*nb*i+nbtx*3
					i++;
					k += *nb;
				}
				if (T1[index1[z]] <= UT[0]) sum[0] -= K[index1[z]]; // compute sum
				else sum[1] += K[index1[z]]; // compute sum
				z++;
			}
			p = sum[0]/SV[index1[z-1]];
			while (k < k0+nbt) { // needed for bootstrap
				P[k+nbtx] -= sum[1]; // k+nbtx = *b+*nb*( i+*nt*j+*nt*(*nx) ) = *b+nbt*j+*nb*i+nbtx
				P[k+nbtx] /= SV[index1[z-1]];
				P[k+nbtx*3] = p; // k+nbtx*3 = *b+*nb*(i+*nt*j+*nt*(*nx)*3) = *b+nbt*j+*nb*i+nbtx*3
				k += *nb;
			}
			while (k > k0) {
				k -= *nb;
				P[k+nbtx] /= P[k0]; // compute and save p12(s,t|X)
				if (P[k+nbtx] < 0) P[k+nbtx] = 0;
				P[k] /= P[k0]; // compute and save p11(s,t|X)
				if (P[k] < 0) P[k] = 0;
				else if (P[k] > 1) P[k] = 1;
				P[k+nbtx*2] = 1-P[k]-P[k+nbtx]; // compute and save p13(s,t|X)
				if (P[k+nbtx*2] < 0) {
					P[k+nbtx] = 1-P[k]; // compute and save p12(s,t|X)
					P[k+nbtx*2] = 0; // save p13(s,t|X)
				}
				P[k+nbtx*3] /= P[k0+nbtx*3]; // compute and save p22(s,t|X)
				if (P[k+nbtx*3] < 0) P[k+nbtx*3] = 0;
				else if (P[k+nbtx*3] > 1) P[k+nbtx*3] = 1;
			}
		}
		free(K); // free memory block
		free(SV); // free memory block
	}
	return;
} // transLIN2I

/*
Author:
	Artur Agostinho Araujo <artur.stat@gmail.com>

Description:
	Computes a conditional transition probability array based
		on the LIN estimator.

Parameters:
	object			an object of class 'LIN2'.
	UT				unique times vector.
	UX				unique covariate vector.
	h				bandwidth parameter.
	window			a string indicating the desired window or kernel.
	methodweights	a string indicating the desired weights method.
	nboot			number of bootstrap samples.

Return value:
	Returns a list where the first element is a
		(nboot)x(nt)x(nx)x4 array of transition probabilities,
		and the second element is the bandwidth value used to
		compute the conditional transition probability estimates.
*/

SEXP TransPROBLIN2(
	SEXP object,
	SEXP UT,
	SEXP UX,
	SEXP h,
	SEXP window,
	SEXP methodweights,
	SEXP nboot)
{
	SEXP data, T1, E1, S, E, X;
	data = VECTOR_ELT(object, 0);
	T1 = VECTOR_ELT(data, 0);
	E1 = VECTOR_ELT(data, 1);
	S = VECTOR_ELT(data, 2);
	E = VECTOR_ELT(data, 3);
	X = VECTOR_ELT(data, 4);
	int len = GET_LENGTH(T1), nt = GET_LENGTH(UT), nx = GET_LENGTH(UX);
	Kfunc kfunc = kchar2ptr(window); // declare and get pointer to function
	Wfunc wfunc = NWWeights; // declare and assign pointer to function
	if (strcmp(CHAR( STRING_ELT(methodweights, 0) ), "LL") == 0) wfunc = LLWeights;
	SEXP dims, P, list;
	PROTECT( dims = allocVector(INTSXP, 4) );
	INTEGER(dims)[0] = *INTEGER(nboot);
	INTEGER(dims)[1] = nt;
	INTEGER(dims)[2] = nx;
	INTEGER(dims)[3] = 4;
	PROTECT( P = allocArray(REALSXP, dims) );
	PROTECT( list = NEW_LIST(2) );
	int b;
	int *index0 = (int*)malloc( len*sizeof(int) ); // allocate memory block
	int *index1 = (int*)malloc( len*sizeof(int) ); // allocate memory block
	double *WORK = (double*)malloc( len*sizeof(double) ); // allocate memory block
	stype SW; // declare stype structure
	SW.type = SINT_PTR; // type is a short int pointer
	SW.ptr.shortinteger = (short int*)malloc( len*sizeof(short int) ); // allocate memory block
	SW.length = len; // hold length of array
	for (b = 0; b < len; b++) SW.ptr.shortinteger[b] = 1; // weights should be equal to 1.0/len, however all weights equal to 1 yield an equivalent result in this case
	b = 0; // b = len, put it back to 0 or a crash might occur
	indx_ii(index0, index1, &len); // initialize indexes
	order_d(REAL(T1), index0, len, FALSE, FALSE, WORK); // get permuation
	order_d(REAL(S), index1, len, FALSE, FALSE, WORK); // get permuation
	transLIN2I(&len, REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), REAL(X), &SW, index0, index1, &nt, REAL(UT), &nx, REAL(UX), REAL(h), kfunc, wfunc, INTEGER(nboot), REAL(P), &b); // compute transition probabilities
	free(index0); // free memory block
	free(index1); // free memory block
	free(WORK); // free memory block
	if (*INTEGER(nboot) > 1) {
		#ifdef _OPENMP
		#pragma omp parallel num_threads(global_num_threads) private(b)
		#endif
		{
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
			#pragma omp for
			#endif
			for (b = 1; b < *INTEGER(nboot); b++) {
				boot_ii(index0, index1, &len); // bootstrap indexes
				order_d(REAL(T1), index0, len, FALSE, FALSE, WORK); // get permuation
				order_d(REAL(S), index1, len, FALSE, FALSE, WORK); // get permuation
				transLIN2I(&len, REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), REAL(X), &SW, index0, index1, &nt, REAL(UT), &nx, REAL(UX), REAL(h), kfunc, wfunc, INTEGER(nboot), REAL(P), &b); // compute transition probabilities
			}
			free(index0); // free memory block
			free(index1); // free memory block
			free(WORK); // free memory block
		}
	}
	free(SW.ptr.shortinteger); // free memory block
	SET_ELEMENT(list, 0, P);
	SET_ELEMENT(list, 1, h);
	UNPROTECT(3);
	return list;
} // TransPROBLIN2
