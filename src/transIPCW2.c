/**************************************/
/*** IPCW2 TRANSITION PROBABILITIES ***/
/**************************************/

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
		p11(s,t|X) = P(Z>t|Z>s|X) = {1-P(Z<=t|X)}/{1-P(Z<=s|X)}
		p12(s,t|X) = P(Z<=t,T>t|Z>s|X) = {P(Z<=t|X)-P(Z<=s|X)-P(s<Z<=t,T<=t|X)}/{1-P(Z<=s|X)}
		p13(s,t|X) = 1-p11(s,t|X)-p12(s,t|X)
		p22(s,t|X) = P(Z<=t,T>t|Z<=s,T>s|X) = {P(Z<=s|X)-P(Z<=s,T<=t|X)}/{P(Z<=s|X)-P(T<=s|X)}

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

static void transIPCW3I(
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
		double sum[4], p[2];
		double *K = (double*)malloc( *len*sizeof(double) ); // allocate memory block
		double *SV = (double*)malloc( *len*sizeof(double) ); // allocate memory block
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for (j = 0; j < *nx; j++) {
			wfunc(X, SW, index0, &UX[j], h, K, kfunc); // compute weights
			wikmsurv(len, T1, E1, K, index0, &e[1], SV); // compute conditional survival probabilities
			for (z = 0, sum[0] = 1; z < e[0]; z++) if (E1[index0[z]] && SV[index0[z]]) sum[0] -= K[index0[z]]/SV[index0[z]]; // compute sum
			for (i = 0, k0 = *b+nbt*j, k = k0, sum[1] = 0; z < e[1]; z++) {
				p[0] = (sum[0]-sum[1])/sum[0];
				p[1] = sum[1]/sum[0];
				while (T1[index0[z]] > UT[i]) {
					if (p[0] < 0) P[k] = 0; // k = *b+*nb*(i+*nt*j) = *b+nbt*j+*nb*i
					else P[k] = p[0];
					P[k+nbtx] = p[1]; // k+nbtx = *b+*nb*( i+*nt*j+*nt*(*nx) ) = *b+nbt*j+*nb*i+nbtx
					i++;
					k += *nb;
				}
				if (E1[index0[z]] && SV[index0[z]]) sum[1] += K[index0[z]]/SV[index0[z]]; // compute sum
			}
			p[0] = (sum[0]-sum[1])/sum[0];
			p[1] = sum[1]/sum[0];
			while (k < k0+nbt) {
				if (p[0] < 0) P[k] = 0; // k = *b+*nb*(i+*nt*j) = *b+nbt*j+*nb*i
				else P[k] = p[0];
				P[k+nbtx] = p[1]; // k+nbtx = *b+*nb*( i+*nt*j+*nt*(*nx) ) = *b+nbt*j+*nb*i+nbtx
				k += *nb;
			}
			wikmsurv(len, S, E, K, index1, &e[3], SV); // compute conditional survival probabilities
			for (z = 0, sum[1] = 0; z < e[2]; z++) if (E[index1[z]] && SV[index1[z]]) sum[1] += K[index1[z]]/SV[index1[z]]; // compute sum
			for (i = 0, k = k0, sum[2] = 0, sum[3] = 0; z < e[3]; z++) {
				p[0] = sum[3]/sum[0];
				p[1] = 1-sum[2]/(1-sum[0]-sum[1]);
				while (S[index1[z]] > UT[i]) {
					P[k+nbtx] -= p[0]; // k+nbtx = *b+*nb*( i+*nt*j+*nt*(*nx) ) = *b+nbt*j+*nb*i+nbtx
					if (P[k+nbtx] < 0) P[k+nbtx] = 0;
					P[k+nbtx*2] = 1-P[k]-P[k+nbtx]; // k+nbtx*2 = *b+*nb*(i+*nt*j+*nt*(*nx)*2) = *b+nbt*j+*nb*i+nbtx*2
					if (P[k+nbtx*2] < 0) {
						P[k+nbtx] = 1-P[k]; // compute and save p12(s,t|X)
						P[k+nbtx*2] = 0; // save p13(s,t|X)
					}
					P[k+nbtx*3] = p[1]; // k+nbtx*3 = *b+*nb*(i+*nt*j+*nt*(*nx)*3) = *b+nbt*j+*nb*i+nbtx*3
					if (P[k+nbtx*3] < 0) P[k+nbtx*3] = 0;
					i++;
					k += *nb;
				}
				if (E[index1[z]] && SV[index1[z]]) {
					if (T1[index1[z]] <= UT[0]) sum[2] += K[index1[z]]/SV[index1[z]]; // compute sum
					else sum[3] += K[index1[z]]/SV[index1[z]]; // compute sum
				}
			}
			p[0] = sum[3]/sum[0];
			p[1] = 1-sum[2]/(1-sum[0]-sum[1]);
			while (k < k0+nbt) {
				P[k+nbtx] -= p[0]; // k+nbtx = *b+*nb*( i+*nt*j+*nt*(*nx) ) = *b+nbt*j+*nb*i+nbtx
				if (P[k+nbtx] < 0) P[k+nbtx] = 0;
				P[k+nbtx*2] = 1-P[k]-P[k+nbtx]; // k+nbtx*2 = *b+*nb*(i+*nt*j+*nt*(*nx)*2) = *b+nbt*j+*nb*i+nbtx*2
				if (P[k+nbtx*2] < 0) {
					P[k+nbtx] = 1-P[k]; // compute and save p12(s,t|X)
					P[k+nbtx*2] = 0; // save p13(s,t|X)
				}
				P[k+nbtx*3] = p[1]; // k+nbtx*3 = *b+*nb*(i+*nt*j+*nt*(*nx)*3) = *b+nbt*j+*nb*i+nbtx*3
				if (P[k+nbtx*3] < 0) P[k+nbtx*3] = 0;
				k += *nb;
			}
		}
		free(K); // free memory block
		free(SV); // free memory block
	}
} // transIPCW3I

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

static void transIPCW4I(
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
	int y, e[4];
	y = 0;
	getIndexI(T1, index0, &UT[0], len, &y, &e[0]); // determine first index
	y = e[0];
	getIndexI(T1, index0, &UT[*nt-1], len, &y, &e[1]); // determine last index
	y = 0;
	getIndexI(S, index1, &UT[0], len, &y, &e[2]); // determine first index
	y = e[2];
	getIndexI(S, index1, &UT[*nt-1], len, &y, &e[3]); // determine last index
	#ifdef _OPENMP
	#pragma omp parallel if(*b < 1) num_threads(global_num_threads) private(y)
	#endif
	{
		register int i, j, z;
		register int64_t k;
		int64_t k0;
		double sum;
		double *K = (double*)malloc( *len*sizeof(double) ); // allocate memory block
		double *SV = (double*)malloc( *len*sizeof(double) ); // allocate memory block
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for (j = 0; j < *nx; j++) {
			wfunc(X, SW, index0, &UX[j], h, K, kfunc); // compute weights
			wikmsurv(len, T1, E1, K, index0, len, SV); // compute conditional survival probabilities
			z = *len-1;
			do {z--;} while (!SV[index0[z]] && z >= e[0]);
			for (sum = 0; z >= e[1]; z--) sum += K[index0[z]]*E1[index0[z]]/SV[index0[z]]; // compute sum
			for (i = *nt-1, k0 = *b+nbt*j, k = k0+nbt-*nb; z >= e[0]; z--) {
				while (T1[index0[z]] <= UT[i]) {
					P[k] = sum; // k = *b+*nb*(i+*nt*j) = *b+nbt*j+*nb*i
					P[k+nbtx] = 0; // k+nbtx = *b+*nb*( i+*nt*j+*nt*(*nx) ) = *b+nbt*j+*nb*i+nbtx
					P[k+nbtx*3] = 0; // k+nbtx*3 = *b+*nb*(i+*nt*j+*nt*(*nx)*3) = *b+nbt*j+*nb*i+nbtx*3
					i--;
					k -= *nb;
				}
				sum += K[index0[z]]*E1[index0[z]]/SV[index0[z]]; // compute sum
			}
			while (k >= k0) {
				P[k] = sum; // k = *b+*nb*(i+*nt*j) = *b+nbt*j+*nb*i
				P[k+nbtx] = 0; // k+nbtx = *b+*nb*( i+*nt*j+*nt*(*nx) ) = *b+nbt*j+*nb*i+nbtx
				P[k+nbtx*3] = 0; // k+nbtx*3 = *b+*nb*(i+*nt*j+*nt*(*nx)*3) = *b+nbt*j+*nb*i+nbtx*3
				k -= *nb;
			}
			wikmsurv(len, S, E, K, index1, len, SV); // compute conditional survival probabilities
			for (z = e[2], y = 0; z < e[3]; z++) {
				while (S[index1[z]] > UT[y]) y++;
				if (E[index1[z]] && SV[index1[z]]) {
					sum = K[index1[z]]/SV[index1[z]];
					if (T1[index1[z]] <= UT[0]) for (k = k0+nbtx*3, i = k+*nb*y; k < i; k += *nb) P[k] += sum; // compute sum
					else for (i = 0, k = k0+nbtx; i < y; i++) P[k+*nb*i] += (T1[index1[z]] <= UT[i])*sum; // compute sum
				}
			}
			for (;z < *len; z++) {
				if (E[index1[z]] && SV[index1[z]]) {
					sum = K[index1[z]]/SV[index1[z]];
					if (T1[index1[z]] <= UT[0]) for (k = k0+nbtx*3, i = k+nbt; k < i; k += *nb) P[k] += sum; // compute sum
					else for (i = 0, k = k0+nbtx; i < *nt; i++) P[k+*nb*i] += (T1[index1[z]] <= UT[i])*sum; // compute sum
				}
			}
			for (k = k0+nbt-*nb; k >= k0; k -= *nb) {
				P[k+nbtx] /= P[k0]; // compute and save p12(s,t|X)
				P[k] /= P[k0]; // compute and save p11(s,t|X)
				P[k+nbtx*2] = 1-P[k]-P[k+nbtx]; // compute and save p13(s,t|X)
				if (P[k+nbtx*2] < 0) {
					P[k+nbtx] = 1-P[k]; // compute and save p12(s,t|X)
					P[k+nbtx*2] = 0; // save p13(s,t|X)
				}
				P[k+nbtx*3] /= P[k0+nbtx*3]; // compute and save p22(s,t|X)
			}
		}
		free(K); // free memory block
		free(SV); // free memory block
	}
} // transIPCW4I

/*
Author:
	Artur Agostinho Araujo <artur.stat@gmail.com>

Description:
	Computes a conditional transition probability array based
		on the inverse probability of censoring estimator.

Parameters:
	object			an object of class 'IPCW2'.
	UT				unique times vector.
	UX				unique covariate vector.
	h				bandwidth parameter.
	window			a string indicating the desired window or kernel.
	methodweights	a string indicating the desired weights method.
	nboot			number of bootstrap samples.
	methodest		an integer indicating the desired method.

Return value:
	Returns a list where the first element is a
		(nboot)x(nt)x(nx)x4 array of transition probabilities,
		and the second element is the bandwidth value used to
		compute the conditional transition probability estimates.
*/

SEXP TransPROBIPCW2(
	SEXP object,
	SEXP UT,
	SEXP UX,
	SEXP h,
	SEXP window,
	SEXP methodweights,
	SEXP nboot,
	SEXP methodest)
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
	void (*func)(CintCP, CdoubleCP, CintCP, CdoubleCP, CintCP, CdoubleCP, CstypeCP, CintCP, CintCP, CintCP, CdoubleCP, CintCP, CdoubleCP, CdoubleCP, Kfunc, Wfunc, CintCP, doubleCP, CintCP);
	switch ( *INTEGER(methodest) ) {
		case 2:
			func = transIPCW4I;
			break;
		default:
			func = transIPCW3I;
	}
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
	func(&len, REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), REAL(X), &SW, index0, index1, &nt, REAL(UT), &nx, REAL(UX), REAL(h), kfunc, wfunc, INTEGER(nboot), REAL(P), &b); // compute transition probabilities
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
				func(&len, REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), REAL(X), &SW, index0, index1, &nt, REAL(UT), &nx, REAL(UX), REAL(h), kfunc, wfunc, INTEGER(nboot), REAL(P), &b); // compute transition probabilities
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
} // TransPROBIPCW2
