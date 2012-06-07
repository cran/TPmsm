/***********************************************/
/*** LOCATION-SCALE TRANSITION PROBABILITIES ***/
/***********************************************/

#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <Rdefines.h>
#include <Rmath.h>
#include "defines.h"
#include "boot.h"
#include "get.h"
#include "sort.h"
#include "splines.h"
#include "unique.h"
#include "window.h"

typedef void (*Kfunc)(CdoubleCP, CintCP, CdoubleCP, CdoubleCP, doubleCP);

#define wsurv_body \
	i = *len-1; \
	n = 0; \
	while (i >= 0) { \
		n += E1[sample0[i]]*W[sample0[i]]; \
		d = E[sample0[i]]*W[sample0[i]]; \
		for (j = i, i--; i >= 0 && T2[sample0[i]] == T2[sample0[i+1]]; i--) { \
			n += E1[sample0[i]]*W[sample0[i]]; \
			d += E[sample0[i]]*W[sample0[i]]; \
		} \
		for (;j > i+1; j--) W[sample0[j]] = 1; \
		if (n != 0) W[sample0[j]] = 1-d/n; \
		else W[sample0[j]] = 1; \
	} \
	for (i = 1; i < *len; i++) { \
		if (T2[sample0[i]] == T2[sample0[i-1]] && W[sample0[i]] != 1) continue; \
		W[sample0[i]] *= W[sample0[i-1]]; \
	} \

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes the kernel bandwidth used by LS estimators
		by cross-validation.

Parameters:
	T1[in]			pointer to T1 first element.
	E1[in]			pointer to E1 first element.
	S[in]			pointer to S first element.
	E[in]			pointer to E first element.
	T2[in]			pointer to T2 first element.
	index[in]		pointer to index first element.
	len[in]			pointer to length of index,
						which must be lower or equal
						than the length of vectors
						T1, E1, S, E and T2.
	h[in]			pointer to h first element.
	nh[in]			pointer to number of bandwidth
						values to test by cross-validation.
	ncv[]			pointer to number of
						cross-validation samples.
	kfunc[in]		pointer to kernel density function.
	b[in]			pointer to row index.
	H[out]			pointer to H.

Return value:
	This function doesn't return a value.

Remarks:
	Vectors T1, E1, S, E and T2 must have the same length.
	Cubic spline interpolation is used.
	If (*b < 1) the code runs in parallel among the available threads.
*/

static void crossValid(
	CdoubleCP T1,
	CintCP E1,
	CdoubleCP S,
	CintCP E,
	CdoubleCP T2,
	CintCP index,
	CintCP len,
	Cdouble h[2],
	CintCP nh,
	CintCP ncv,
	Kfunc kfunc,
	CintCP b,
	doubleCP H)
{
	if (h[1] == h[0]) {
		*H = h[0];
		return;
	}
	register int x, y, z, i, j;
	const int method = 1;
	int u0;
	double h0, aux[3], cv0, cv1, cv2, n, d;
	h0 = (h[1]-h[0])/(*nh-1);
	cv0 = R_PosInf;
	#ifdef _OPENMP
	#pragma omp parallel if(*b < 1) private(x, y, z, i, j, u0, aux, cv1, cv2, n, d)
	#endif
	{
		double *h1 = (double*)malloc( sizeof(double) ); // allocate memory block
		int *sample0 = (int*)malloc( *len*sizeof(int) ); // allocate memory block
		int *sample1 = (int*)malloc( *len*sizeof(int) ); // allocate memory block
		double *W = (double*)malloc( *len*sizeof(double) ); // allocate memory block
		double *MX = (double*)malloc( *len*sizeof(double) ); // allocate memory block
		int *unique0 = (int*)malloc( *len*sizeof(int) ); // allocate memory block
		double *a = (double*)malloc( *len*3*sizeof(double) ); // allocate memory block
		if (*b < 1) {
			unsigned int iseed = (unsigned int)time(NULL) ^ (unsigned int)omp_get_thread_num(); // save per thread seed
			srand(iseed); // set seed
		}
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for (x = 0; x < *nh; x++) {
			*h1 = h[0]+h0*x;
			for (cv1 = 0, y = 0; y < *ncv; y++) {
				boot_i(sample0, len); // simulate first random sample
				boot_i(sample1, len); // simulate second random sample
				for (z = 0; z < *len; z++) {
					sample0[z] = index[sample0[z]];
					sample1[z] = index[sample1[z]];
				}
				order_d(T1, sample0, *len, FALSE, FALSE, MX); // get permuation
				uniqueI(T1, sample0, len, unique0, &u0); // compute unique index
				order_d(T2, sample0, *len, FALSE, FALSE, MX); // get permuation
				for (z = 0; z < u0; z++) {
					kfunc(T1, len, &T1[unique0[z]], h1, W); // compute weights
					wsurv_body // compute survival probabilities vector
					MX[unique0[z]] = (1-W[sample0[0]])*T2[sample0[0]]; // initialize mean
					for (i = 1; i < *len; i++) {
						MX[unique0[z]] += (W[sample0[i-1]]-W[sample0[i]])*T2[sample0[i]]; // compute mean
					}
					if (W[sample0[*len-1]] != 1) MX[unique0[z]] /= 1-W[sample0[*len-1]]; // normalize mean
				}
				spline_coefI(&method, T1, MX, unique0, &u0, a, a+*len, a+*len*2); // compute spline coefficients
				spline_evalI(&method, T1, MX, unique0, &u0, a, a+*len, a+*len*2, T1, W, sample1, len); // interpolate mean from unique0 to sample1
				order_di(S, E, sample1, *len, FALSE, FALSE, TRUE, MX, unique0); // get permuation
				for (aux[0] = 1, d = 0, cv2 = 0, z = 0; z < *len; z++) {
					aux[2] = (double)E[sample1[z]]/(*len-z); // compute needed factor
					aux[1] = 1-aux[2]; // factor needed for the computation
					aux[2] *= aux[0]; // compute and save weight
					aux[0] *= aux[1]; // compute and save factor needed for next iteration
					d += aux[2]; // sum weights
					cv2 += aux[2]*R_pow_di(T2[sample1[z]]-W[sample1[z]]*(W[sample1[z]] > 0), 2);
				}
				cv2 /= d; // compute cross-validation error term
				cv1 += cv2 / *ncv; // compute cross-validation error
			}
			#ifdef _OPENMP
			#pragma omp critical
			#endif
			{
				if (cv1 < cv0) {
					cv0 = cv1;
					*H = *h1;
				}
			}
		}
		free(h1); // free memory block
		free(sample0); // free memory block
		free(sample1); // free memory block
		free(W); // free memory block
		free(MX); // free memory block
		free(unique0); // free memory block
		free(a); // free memory block
	}
	return;
} // crossValid

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes the mean and variance vectors respectively
		labeled MX and SX.

Parameters:
	T1[in]			pointer to T1 first element.
	E1[in]			pointer to E1 first element.
	T2[in]			pointer to T2 first element.
	E[in]			pointer to E first element.
	index[in]		pointer to index first element.
	len[in]			pointer to length of index,
						which must be lower or equal
						than the length of vectors
						T1, E1, T2, E, MX and SX.
	H[in]			pointer to H.
	kfunc[in]		pointer to kernel density function.
	MX[out]			pointer to MX vector.
	SX[out]			pointer to SX vector.
	b[in]			pointer to row index.

Return value:
	This function doesn't return a value.

Remarks:
	Vector index must indicate the permutation of vector T2
		sorted by ascending order.
	Vectors T1, E1, T2, E, MX and SX must have the same length.
	If (*b < 1) the code runs in parallel among the available threads.
*/

static void LSmeasuresI(
	CdoubleCP T1,
	CintCP E1,
	CdoubleCP T2,
	CintCP E,
	CintCP index,
	CintCP len,
	CdoubleCP H,
	Kfunc kfunc,
	doubleCP MX,
	doubleCP SX,
	CintCP b)
{
	#ifdef _OPENMP
	#pragma omp parallel if(*b < 1)
	#endif
	{
		register int z, i, j;
		double n, d;
		double *W = (double*)malloc( *len*sizeof(double) ); // allocate memory block
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for (z = 0; z < *len; z++) {
			kfunc(T1, len, &T1[index[z]], H, W); // compute weights
			#define sample0 index
			wsurv_body // compute survival probabilities vector
			#undef sample0
			d = 1-W[index[0]];
			MX[index[z]] = d*T2[index[0]]; // initialize mean
			SX[index[z]] = d*R_pow_di(T2[index[0]], 2); // initialize variance
			for (i = 1; i < *len; i++) {
				d = W[index[i-1]]-W[index[i]]; // compute survival probability jump
				MX[index[z]] += d*T2[index[i]]; // compute mean
				SX[index[z]] += d*R_pow_di(T2[index[i]], 2); // compute variance
			}
			d = 1;
			if (W[index[*len-1]] != 1) d /= 1-W[index[*len-1]];
			MX[index[z]] *= d; // normalize mean
			SX[index[z]] *= d; // normalize variance
			SX[index[z]] -= R_pow_di(MX[index[z]], 2); // compute variance
			if (SX[index[z]] < 1e-10) SX[index[z]] = 1e-10; // variance can't be negative
			SX[index[z]] = sqrt(SX[index[z]]); // compute standard deviation
		}
		free(W); // free memory block
	}
	return;
} // LSmeasuresI

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes the transition probabilities:
		p11(s,t) = P(Z>t|Z>s) = P(Z>t)/P(Z>s)
		p12(s,t) = P(Z<=t,T>t|Z>s) = P(s<Z<=t,T>t)/P(Z>s)
		p22(s,t) = P(Z<=t,T>t|Z<=s,T>s) = P(Z<=s,T>t)/P(Z<=s,T>s)

Parameters:
	len[in]			pointer to length of T1, E1, T2 and E.
	T1[in]			pointer to T1 first element.
	E1[in]			pointer to E1 first element.
	T2[in]			pointer to T2 first element.
	E[in]			pointer to E first element.
	index0[in]		pointer to index0 first element.
	index1[inout]	pointer to index1 first element.
	nt[in]			pointer to length of UT and number of rows of P.
	UT[in]			pointer to unique times vector.
	nb[in]			pointer to number of rows of P.
	P[out]			pointer to a (nb)x(nt)x3 probability array.
	b[in]			pointer to row index.
	kfunc[in]		pointer to kernel density function.
	H[in]			pointer to H.

Return value:
	This function doesn't return a value.

Remarks:
	Vector index0 must indicate the permutation of vector T1
		sorted by ascending order.
	Vector index1 must indicate the permutation of vector T2
		sorted by ascending order.
	Vectors T1, E1, T2 and E must have the same length.
*/

static void transLSI(
	CintCP len,
	Cdouble T1[*len],
	Cint E1[*len],
	Cdouble T2[*len],
	Cint E[*len],
	Cint index0[*len],
	int index1[*len],
	CintCP nt,
	Cdouble UT[*nt],
	CintCP nb,
	double P[*nb*(*nt)*3],
	CintCP b,
	Kfunc kfunc,
	CdoubleCP H)
{
	register int i;
	int j, k, x, n, d;
	double aux[2], p;
	double *MX = (double*)malloc( *len*sizeof(double) ); // allocate memory block
	double *SX = (double*)malloc( *len*sizeof(double) ); // allocate memory block
	LSmeasuresI(T1, E1, T2, E, index1, len, H, kfunc, MX, SX, b); // compute location-scale measures
	double *EX = (double*)malloc( *len*sizeof(double) ); // allocate memory block
	double *SURV = (double*)malloc( *len*sizeof(double) ); // allocate memory block
	double *W = (double*)malloc( *len*sizeof(double) ); // allocate memory block
	for (i = 0; i < *len; i++) EX[index1[i]] = (T2[index1[i]]-MX[index1[i]])/SX[index1[i]];
	order_d(EX, index1, *len, FALSE, FALSE, W); // get permuation
	i = 0;
	p = 1;
	while (i < *len) { // loop through the sample until last index is reached
		n = *len-i; // count the living
		d = E[index1[i]]; // initialize dead count
		for (i++; i < *len && EX[index1[i]] == EX[index1[i-1]]; i++) { // loop until time changes or last index is reached
			d += E[index1[i]]; // count the dead
		}
		p *= 1-(double)d/n;
		for (j = *len-n; j < i; j++) SURV[index1[j]] = p;
	}
	j = 0;
	getIndexI(T1, index0, &UT[0], len, &j, &n); // determine first index
	for (aux[0] = 1, p = 1; j < n; j++) { // loop through the sample until last index is reached
		W[index0[j]] = (double)E1[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-W[index0[j]]; // factor needed for the computation
		W[index0[j]] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		p -= W[index0[j]]; // compute probability
	}
	getIndexI(T1, index0, &UT[*nt-1], len, &j, &d); // determine last index
	for (i = 0; j < d; j++) { // loop through the sample until last index is reached
		while (T1[index0[j]] > UT[i]) {
			P[*b+*nb*i] = p; // save transition probability
			i++;
		}
		W[index0[j]] = (double)E1[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-W[index0[j]]; // factor needed for the computation
		W[index0[j]] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		p -= W[index0[j]]; // compute probability
	}
	for (; i < *nt; i++) P[*b+*nb*i] = p;
	x = n;
	#ifdef _OPENMP
	#pragma omp parallel for if(*b < 1) private(i, j, k, aux, d) firstprivate(x) ordered
	#endif
	for (i = 0; i < *nt; i++) {
		#ifdef _OPENMP
		#pragma omp ordered
		#endif
		{
			for (j = 0, P[*b+*nb*(i+*nt*2)] = 0; j < n; j++) {
				k = 0;
				aux[0] = (UT[i]-T1[index0[j]]-MX[index0[j]])/SX[index0[j]];
				aux[1] = 1;
				getOrdinateI(EX, SURV, index1, len, &k, &aux[0], &aux[1]); // get survival probability
				P[*b+*nb*(i+*nt*2)] += aux[1]*W[index0[j]];
			}
			getIndexI(T1, index0, &UT[i], len, &x, &d); // determine last index
			for (P[*b+*nb*(i+*nt)] = 0; j < d; j++) {
				k = 0;
				aux[0] = (UT[i]-T1[index0[j]]-MX[index0[j]])/SX[index0[j]];
				aux[1] = 1;
				getOrdinateI(EX, SURV, index1, len, &k, &aux[0], &aux[1]); // get survival probability
				P[*b+*nb*(i+*nt)] += aux[1]*W[index0[j]];
			}
			x = d;
		}
	}
	free(MX); // free memory block
	free(SX); // free memory block
	free(EX); // free memory block
	free(SURV); // free memory block
	free(W); // free memory block
	#ifdef _OPENMP
	#pragma omp parallel for if(*b < 1) private(i) ordered
	#endif
	for (i = *nt-1; i >= 0; i--) {
		#ifdef _OPENMP
		#pragma omp ordered
		#endif
		{
			P[*b+*nb*(i+*nt)] /= P[*b];
			P[*b+*nb*i] /= P[*b];
			P[*b+*nb*(i+*nt*2)] /= P[*b+*nb*(*nt*2)];
		}
	}
	return;
} // transLSI

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes a transition probability array
		based on the Location-Scale estimator.

Parameters:
	object			an object of class 'LS'.
	UT				unique times vector.
	h				a vector of bandwidth values with length two.
	nh				number of bandwidth values
						to test by cross-validation.
	ncv				number of cross-validation samples.
	window			a string indicating the desired window or kernel.
	nboot			number of bootstrap samples.

Return value:
	Returns a list where the first element is a
		(nboot)x(nt)x3 array of transition probabilities,
		and the second element is the bandwidth value used
		to compute the transition probability estimates.
*/

SEXP TransPROBLS(
	SEXP object,
	SEXP UT,
	SEXP h,
	SEXP nh,
	SEXP ncv,
	SEXP window,
	SEXP nboot,
	SEXP bootcv)
{
	SEXP data, T1, E1, S, E;
	data = VECTOR_ELT(object, 0);
	T1 = VECTOR_ELT(data, 0);
	E1 = VECTOR_ELT(data, 1);
	S = VECTOR_ELT(data, 2);
	E = VECTOR_ELT(data, 3);
	int len = GET_LENGTH(T1), nt = GET_LENGTH(UT), b;
	double *T2 = (double*)malloc( len*sizeof(double) ); // allocate memory block
	for (b = 0; b < len; b++) T2[b] = REAL(S)[b]-REAL(T1)[b];
	Kfunc kfunc; // declare pointer to function
	kfunc = kchar2ptr(window); // get pointer to function
	SEXP P, H, list;
	PROTECT( P = alloc3DArray(REALSXP, *INTEGER(nboot), nt, 3) );
	PROTECT( H = NEW_NUMERIC(1) );
	PROTECT( list = NEW_LIST(2) );
	int *index0 = (int*)malloc( len*sizeof(int) ); // allocate memory block
	int *index1 = (int*)malloc( len*sizeof(int) ); // allocate memory block
	double *WORK0 = (double*)malloc( len*sizeof(double) ); // allocate memory block
	int *WORK1 = (int*)malloc( len*sizeof(int) ); // allocate memory block
	double HC;
	b = 0;
	indx_ii(index0, index1, &len); // initialize indexes
	crossValid(REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), T2, index0, &len, REAL(h), INTEGER(nh), INTEGER(ncv), kfunc, &b, &HC); // compute bandwidth
	order_di(REAL(T1), INTEGER(E1), index0, len, FALSE, FALSE, TRUE, WORK0, WORK1); // get permuation
	order_d(T2, index1, len, FALSE, FALSE, WORK0); // get permuation
	transLSI(&len, REAL(T1), INTEGER(E1), T2, INTEGER(E), index0, index1, &nt, REAL(UT), INTEGER(nboot), REAL(P), &b, kfunc, &HC); // compute transition probabilities
	free(index0); // free memory block
	free(index1); // free memory block
	free(WORK0); // free memory block
	free(WORK1); // free memory block
	*REAL(H) = HC;
	if (*INTEGER(nboot) > 1) {
		#ifdef _OPENMP
		#pragma omp parallel private(b) firstprivate(HC)
		#endif
		{
			int *index0 = (int*)malloc( len*sizeof(int) ); // allocate memory block
			int *index1 = (int*)malloc( len*sizeof(int) ); // allocate memory block
			double *WORK0 = (double*)malloc( len*sizeof(double) ); // allocate memory block
			int *WORK1 = (int*)malloc( len*sizeof(int) ); // allocate memory block
			unsigned int iseed = (unsigned int)time(NULL) ^ (unsigned int)omp_get_thread_num(); // save per thread seed
			srand(iseed); // set seed
			#ifdef _OPENMP
			#pragma omp for
			#endif
			for (b = 1; b < *INTEGER(nboot); b++) {
				boot_ii(index0, index1, &len); // bootstrap indexes
				if ( *LOGICAL(bootcv) ) crossValid(REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), T2, index0, &len, REAL(h), INTEGER(nh), INTEGER(ncv), kfunc, &b, &HC); // compute bandwidth
				order_di(REAL(T1), INTEGER(E1), index0, len, FALSE, FALSE, TRUE, WORK0, WORK1); // get permuation
				order_d(T2, index1, len, FALSE, FALSE, WORK0); // get permuation
				transLSI(&len, REAL(T1), INTEGER(E1), T2, INTEGER(E), index0, index1, &nt, REAL(UT), INTEGER(nboot), REAL(P), &b, kfunc, &HC); // compute transition probabilities
			}
			free(index0); // free memory block
			free(index1); // free memory block
			free(WORK0); // free memory block
			free(WORK1); // free memory block
		}
	}
	free(T2); // free memory block
	SET_ELEMENT(list, 0, P);
	SET_ELEMENT(list, 1, H);
	UNPROTECT(3);
	return list;
} // TransPROBLS
