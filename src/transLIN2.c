/*************************************/
/*** LIN2 TRANSITION PROBABILITIES ***/
/*************************************/

#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <Rdefines.h>
#include "defines.h"
#include "boot.h"
#include "get.h"
#include "sort.h"
#include "window.h"

typedef void (*Kfunc)(CdoubleCP, CintCP, CdoubleCP, CdoubleCP, doubleCP);
typedef void (*Wfunc)(CdoubleCP, CintCP, CdoubleCP, CdoubleCP, doubleCP, Kfunc kfunc);

#define wisurv_body \
	for (i = *len-1, n = 0; i >= e; i--) { \
		n += W[index0[i]]; \
	} \
	while (i >= 0) { \
		n += W[index0[i]]; \
		r = E1[index0[i]]*W[index0[i]]; \
		d = (1-E1[index0[i]])*W[index0[i]]; \
		for (k = i, i--; i >= 0 && T1[index0[i]] == T1[index0[i+1]]; i--) { \
			n += W[index0[i]]; \
			r += E1[index0[i]]*W[index0[i]]; \
			d += (1-E1[index0[i]])*W[index0[i]]; \
		} \
		for (;k > i+1; k--) SURV[index0[k]] = 1; \
		if (n-r != 0) SURV[index0[k]] = 1-d/(n-r); \
		else SURV[index0[k]] = 1; \
	} \
	for (i = 1; i < e; i++) { \
		if (T1[index0[i]] == T1[index0[i-1]] && SURV[index0[i]] != 1) continue; \
		SURV[index0[i]] *= SURV[index0[i-1]]; \
	} \

/*
Author:
	Artur Agostinho Ara�jo <b5498@math.uminho.pt>

Description:
	Computes the conditional transition probabilities:
		p11(s,t|X) = P(Z>t|Z>s|X) = P(Z>t|X)/P(Z>s|X)
		p12(s,t|X) = P(Z<=t,T>t|Z>s|X) = P(s<Z<=t,T>t|X)/P(Z>s|X)
		p22(s,t|X) = P(Z<=t,T>t|Z<=s,T>s|X) = P(Z<=s,T>t|X)/P(Z<=s,T>s|X)

Parameters:
	len[in]			pointer to length of T1, E1, S, E and X.
	T1[in]			pointer to T1 first element.
	E1[in]			pointer to E1 first element.
	S[in]			pointer to S first element.
	E[in]			pointer to E first element.
	X[in]			pointer to X first element.
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
	P[out]			pointer to a (nb)x(nt)x(nx)x3 probability array.
	b[in]			pointer to row index.

Return value:
	This function doesn't return a value.

Remarks:
	Vector index0 must indicate the permutation of vector T1
		sorted by ascending order.
	Vector index1 must indicate the permutation of vector S
		sorted by ascending order.
	Vectors T1, E1, S, E and X must have the same length.
*/

static void transLIN2I(
	CintCP len,
	Cdouble T1[*len],
	Cint E1[*len],
	Cdouble S[*len],
	Cint E[*len],
	Cdouble X[*len],
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
	double P[*nb*(*nt)*(*nx)*3],
	CintCP b)
{
	const int nrc = (*nt)*(*nx);
	int k, e[4];
	k = 0;
	getIndexI(T1, index0, &UT[0], len, &k, &e[0]); // determine first index
	k = e[0];
	getIndexI(T1, index0, &UT[*nt-1], len, &k, &e[1]); // determine last index
	k = 0;
	getIndexI(S, index1, &UT[0], len, &k, &e[2]); // determine first index
	k = e[2];
	getIndexI(S, index1, &UT[*nt-1], len, &k, &e[3]); // determine last index
	#ifdef _OPENMP
	#pragma omp parallel if(*b < 1) private(k)
	#endif
	{
		register int i, j;
		double p, n, r, d, sum[2];
		double *W = (double*)malloc( *len*sizeof(double) ); // allocate memory block
		double *SURV = (double*)malloc( *len*sizeof(double) ); // allocate memory block
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for (j = 0; j < *nx; j++) {
			wfunc(X, len, &UX[j], h, W, kfunc); // compute weights
			#define e e[1]
			wisurv_body // compute survival probabilities
			#undef e
			for (k = 0, sum[0] = 0; k < e[0]; k++) sum[0] += W[index0[k]];
			i = 0; sum[1] = 0;
			if (k == 0) {
				while (T1[index0[0]] > UT[i]) {
					P[*b+*nb*(i+*nt*j)] = n;
					P[*b+*nb*(i+*nt*j+nrc)] = 0;
					i++;
				}
				sum[1] += W[index0[0]];
				k++;
			}
			while (k < e[1]) {
				p = (n-sum[0]-sum[1])/SURV[index0[k-1]];
				while (T1[index0[k]] > UT[i]) {
					P[*b+*nb*(i+*nt*j)] = p;
					P[*b+*nb*(i+*nt*j+nrc)] = sum[1];
					i++;
				}
				sum[1] += W[index0[k]];
				k++;
			}
			p = (n-sum[0]-sum[1])/SURV[index0[k-1]];
			for (; i < *nt; i++) { // needed for bootstrap
				P[*b+*nb*(i+*nt*j)] = p;
				P[*b+*nb*(i+*nt*j+nrc)] = sum[1];
			}
			#define e e[3]
			#define E1 E
			#define index0 index1
			#define T1 S
			wisurv_body // compute survival probabilities
			#undef e
			#undef E1
			#undef index0
			#undef T1
			for (k = 0, sum[1] = 0; k < e[2]; k++) sum[1] += W[index1[k]];
			i = 0; sum[0] -= sum[1]; sum[1] = 0;
			if (k == 0) {
				while (S[index1[0]] > UT[i]) {
					P[*b+*nb*(i+*nt*j+nrc*2)] = sum[0];
					i++;
				}
				if (T1[index1[0]] <= UT[0]) sum[0] -= W[index1[0]]; // compute sum
				else sum[1] += W[index1[0]]; // compute sum
				k++;
			}
			while (k < e[3]) {
				p = sum[0]/SURV[index1[k-1]];
				while (S[index1[k]] > UT[i]) {
					P[*b+*nb*(i+*nt*j+nrc)] -= sum[1];
					P[*b+*nb*(i+*nt*j+nrc)] /= SURV[index1[k-1]];
					P[*b+*nb*(i+*nt*j+nrc*2)] = p;
					i++;
				}
				if (T1[index1[k]] <= UT[0]) sum[0] -= W[index1[k]]; // compute sum
				else sum[1] += W[index1[k]]; // compute sum
				k++;
			}
			p = sum[0]/SURV[index1[k-1]];
			for (; i < *nt; i++) { // needed for bootstrap
				P[*b+*nb*(i+*nt*j+nrc)] -= sum[1];
				P[*b+*nb*(i+*nt*j+nrc)] /= SURV[index1[k-1]];
				P[*b+*nb*(i+*nt*j+nrc*2)] = p;
			}
			for (i = *nt-1; i >= 0; i--) {
				P[*b+*nb*(i+*nt*j+nrc)] /= P[*b+*nb*(*nt*j)];
				P[*b+*nb*(i+*nt*j)] /= P[*b+*nb*(*nt*j)];
				P[*b+*nb*(i+*nt*j+nrc*2)] /= P[*b+*nb*(*nt*j+nrc*2)];
			}
		}
		free(W); // free memory block
		free(SURV); // free memory block
	}
	return;
} // transLIN2I

/*
Author:
	Artur Agostinho Ara�jo <b5498@math.uminho.pt>

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
		(nboot)x(nt)x(nx)x3 array of transition probabilities,
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
	Kfunc kfunc; // declare pointer to function
	kfunc = kchar2ptr(window); // get pointer to function
	Wfunc wfunc; // declare pointer to function
	if (strcmp(CHAR( STRING_ELT(methodweights, 0) ), "LL") == 0) wfunc = LLWeights;
	else wfunc = NWWeights;
	SEXP dims, P, list;
	PROTECT( dims = allocVector(INTSXP, 4) );
	INTEGER(dims)[0] = *INTEGER(nboot);
	INTEGER(dims)[1] = nt;
	INTEGER(dims)[2] = nx;
	INTEGER(dims)[3] = 3;
	PROTECT( P = allocArray(REALSXP, dims) );
	PROTECT( list = NEW_LIST(2) );
	int b = 0;
	int *index0 = (int*)malloc( len*sizeof(int) ); // allocate memory block
	int *index1 = (int*)malloc( len*sizeof(int) ); // allocate memory block
	double *WORK = (double*)malloc( len*sizeof(double) ); // allocate memory block
	indx_ii(index0, index1, &len); // initialize indexes
	order_d(REAL(T1), index0, len, FALSE, FALSE, WORK); // get permuation
	order_d(REAL(S), index1, len, FALSE, FALSE, WORK); // get permuation
	transLIN2I(&len, REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), REAL(X), index0, index1, &nt, REAL(UT), &nx, REAL(UX), REAL(h), kfunc, wfunc, INTEGER(nboot), REAL(P), &b); // compute transition probabilities
	free(index0); // free memory block
	free(index1); // free memory block
	free(WORK); // free memory block
	if (*INTEGER(nboot) > 1) {
		#ifdef _OPENMP
		#pragma omp parallel private(b)
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
				transLIN2I(&len, REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), REAL(X), index0, index1, &nt, REAL(UT), &nx, REAL(UX), REAL(h), kfunc, wfunc, INTEGER(nboot), REAL(P), &b); // compute transition probabilities
			}
			free(index0); // free memory block
			free(index1); // free memory block
			free(WORK); // free memory block
		}
	}
	SET_ELEMENT(list, 0, P);
	SET_ELEMENT(list, 1, h);
	UNPROTECT(3);
	return list;
} // TransPROBLIN2
