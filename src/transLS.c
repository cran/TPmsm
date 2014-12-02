/***********************************************/
/*** LOCATION-SCALE TRANSITION PROBABILITIES ***/
/***********************************************/

#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdlib.h>
#include <Rdefines.h>
#include <Rmath.h>
#include "defines.h"
#include "get.h"
#include "kmsurv.h"
#include "RngStream.h"
#include "RngArray.h"
#include "RngBoot.h"
#include "rthreads.h"
#include "sort.h"
#include "splines.h"
#include "unique.h"
#include "wkmsurv.h"
#include "wtypefunc.h"
#include "window.h"

typedef struct {
	int *sample0, *sample1, *unique0;
	double *K, *MX;
	double *a, *b, *c;
	double *E0B, *E1B;
} transLSW;

#define weights \
	order_di(T1, E1, sample0, *len, FALSE, FALSE, TRUE, a, unique0); /* get permuation */ \
	for (aux[0] = 1, i = 0; i < *len; i++) { /* loop through the sample until last index is reached */ \
		W[sample0[i]] = (double)E1[sample0[i]]/(*len-i); /* compute needed factor */ \
		aux[1] = 1-W[sample0[i]]; /* factor needed for the computation */ \
		W[sample0[i]] *= aux[0]; /* compute and save weight */ \
		aux[0] *= aux[1]; /* compute and save factor needed for next iteration */ \
	} \

#define mean \
	order_d(T2, sample0, *len, FALSE, FALSE, a); /* get permuation */ \
	for (i = 0; i < u0; i++) { \
		kfunc(T1, &SW, sample0, &T1[unique0[i]], &h1, K); /* compute weights */ \
		wkmsurv(len, T2, E, K, sample0, len, K); /* compute conditional survival probabilities vector */ \
		MX[unique0[i]] = (1-K[sample0[0]])*T2[sample0[0]]; /* initialize mean */ \
		for (j = 1; j < *len; j++) { \
			MX[unique0[i]] += (K[sample0[j-1]]-K[sample0[j]])*T2[sample0[j]]; /* compute mean */ \
		} \
		if (K[sample0[*len-1]] != 1) MX[unique0[i]] /= 1-K[sample0[*len-1]]; /* normalize mean */ \
	} \

#define cverror	\
	order_di(S, E, sample1, *len, FALSE, FALSE, TRUE, a, unique0); /* get permuation */ \
	for (aux[0] = 1, sum = 0, cv2 = 0, i = 0; i < *len; i++) { \
		aux[2] = (double)E[sample1[i]]/(*len-i); /* compute needed factor */ \
		aux[1] = 1-aux[2]; /* factor needed for the computation */ \
		aux[2] *= aux[0]; /* compute and save weight */ \
		aux[0] *= aux[1]; /* compute and save factor needed for next iteration */ \
		sum += aux[2]; /* sum weights */ \
		cv2 += aux[2]*R_pow_di(T2[sample1[i]]-K[sample1[i]]*(K[sample1[i]] > 0), 2); \
	} \
	cv2 /= sum; /* compute cross-validation error term */ \

/*
Author:
	Artur Araujo <artur.stat@gmail.com>

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
	ncv[in]			pointer to number of
					cross-validation samples.
	cvfull[in]		pointer to integer indicating if full
					cross-validation is to be done.
	kfunc[in]		pointer to kernel density function.
	H[out]			pointer to H first element.
	t[in]			pointer to thread number.
	WORK[out]		pointer to array of transLSW structures.

Return value:
	This function doesn't return a value.

Remarks:
	Vectors T1, E1, S, E and T2 must have the same length.
	Cubic spline interpolation is used.
*/

static void crossValid(
	CdoubleCP T1,
	intCP E1,
	CdoubleCP S,
	CintCP E,
	CdoubleCP T2,
	CintCP index,
	CintCP len,
	Cdouble h[4],
	CintCP nh,
	CintCP ncv,
	CintCP cvfull,
	Kfunc kfunc,
	double H[2],
	CintCP t,
	transLSW *WORK)
{
	if (h[1] == h[0] && h[3] == h[2]) {
		H[0] = h[0];
		H[1] = h[2];
		return;
	}
	register int x, y, i, j;
	const int method = 1;
	unsigned long iseed[6];
	int u0, tid = *t;
	double h0, h1, aux[3], cv0, cv1, cv2, sum;
	stype SW;
	SW.type = INT_PTR; // type is an int pointer
	SW.ptr.integer = E1; // hold E1 pointer in ptr union
	SW.length = *len; // hold length of array
	#ifdef _OPENMP
	#pragma omp parallel if( !omp_in_parallel() ) num_threads(global_num_threads) firstprivate(tid) private(x, y, i, j, iseed, u0, h1, aux, cv1, cv2, sum)
	#endif
	{
		#ifdef _OPENMP
		if (omp_get_num_threads() != 1) tid = omp_get_thread_num(); // use the correct thread number
		#endif
		int *sample0 = WORK[tid].sample0;
		int *sample1 = WORK[tid].sample1;
		double *K = WORK[tid].K;
		double *MX = WORK[tid].MX;
		int *unique0 = WORK[tid].unique0;
		double *a = WORK[tid].a;
		double *b = WORK[tid].b;
		double *c = WORK[tid].c;
		RngStream_GetState(RngArray[tid], iseed); // save per thread seed
		if (h[1] == h[0]) {
			#ifdef _OPENMP
			#pragma omp single
			#endif
			{
				H[0] = h[0];
			}
		} else {
			#ifdef _OPENMP
			#pragma omp single
			#endif
			{
				h0 = (h[1]-h[0])/(*nh-1);
				cv0 = R_PosInf;
			}
			#ifdef _OPENMP
			#pragma omp for
			#endif
			for (x = 0; x < *nh; x++) {
				h1 = h[0]+h0*x;
				for (cv1 = 0, y = 0; y < *ncv; y++) {
					boot_i(RngArray[tid], len, sample0); // simulate first random sample
					boot_i(RngArray[tid], len, sample1); // simulate second random sample
					for (i = 0; i < *len; i++) {
						sample0[i] = index[sample0[i]];
						sample1[i] = index[sample1[i]];
					}
					order_d(T1, sample0, *len, FALSE, FALSE, a); // get permuation
					uniqueI(T1, sample0, len, unique0, &u0); // compute unique index
					mean // compute mean vector
					spline_coefI(&method, T1, MX, unique0, &u0, a, b, c); // compute spline coefficients
					spline_evalI(&method, T1, MX, unique0, &u0, a, b, c, T1, K, sample1, len); // interpolate mean from unique0 to sample1
					cverror // compute cross-validation error term
					cv1 += cv2 / *ncv; // compute cross-validation error
				}
				#ifdef _OPENMP
				#pragma omp critical
				#endif
				{
					if (cv1 < cv0) {
						cv0 = cv1;
						H[0] = h1;
					}
				}
			}
		}
		if (*cvfull) {
			if (h[3] == h[2]) {
				#ifdef _OPENMP
				#pragma omp single
				#endif
				{
					H[1] = h[2];
				}
			} else {
				RngStream_SetSeed(RngArray[tid], iseed); // restore seed
				double *E0B = WORK[tid].E0B;
				double *E1B = WORK[tid].E1B;
				#ifdef _OPENMP
				#pragma omp single
				#endif
				{
					h0 = (h[3]-h[2])/(*nh-1);
					cv0 = R_PosInf;
				}
				#ifdef _OPENMP
				#pragma omp for
				#endif
				for (x = 0; x < *nh; x++) {
					h1 = h[2]+h0*x;
					for (cv1 = 0, y = 0; y < *ncv; y++) {
						boot_i(RngArray[tid], len, sample0); // simulate first random sample
						boot_i(RngArray[tid], len, sample1); // simulate second random sample
						for (i = 0; i < *len; i++) {
							sample0[i] = index[sample0[i]];
							sample1[i] = index[sample1[i]];
						}
						order_d(T1, sample0, *len, FALSE, FALSE, a); // get permuation
						uniqueI(T1, sample0, len, unique0, &u0); // compute unique index
						#define h1 H[0]
						mean // compute mean vector
						#undef h1
						spline_coefI(&method, T1, MX, unique0, &u0, a, b, c); // compute spline coefficients
						spline_evalI(&method, T1, MX, unique0, &u0, a, b, c, T1, MX, sample0, len); // interpolate mean from unique0 to sample0
						spline_evalI(&method, T1, MX, unique0, &u0, a, b, c, T1, K, sample1, len); // interpolate mean from unique0 to sample1
						for (i = 0; i < *len; i++) {
							E0B[sample0[i]] = T2[sample0[i]]-MX[sample0[i]]; // compute error
							E1B[sample1[i]] = T2[sample1[i]]-K[sample1[i]]; // compute error
						}
						#define T2 E0B
						mean // compute mean vector
						#undef T2
						spline_coefI(&method, T1, MX, unique0, &u0, a, b, c); // compute spline coefficients
						spline_evalI(&method, T1, MX, unique0, &u0, a, b, c, T1, K, sample1, len); // interpolate mean from unique0 to sample1
						for (i = 0; i < *len; i++) MX[sample1[i]] = T1[sample1[i]]+E1B[sample1[i]]; // compute total time
						#define S MX
						#define T2 E1B
						cverror // compute cross-validation error term
						#undef S
						#undef T2
						cv1 += cv2 / *ncv; // compute cross-validation error
					}
					#ifdef _OPENMP
					#pragma omp critical
					#endif
					{
						if (cv1 < cv0) {
							cv0 = cv1;
							H[1] = h1;
						}
					}
				}
			}
		} else {
			#ifdef _OPENMP
			#pragma omp single
			#endif
			{
				H[1] = H[0];
			}
		}
	}
	return;
} // crossValid

/*
Author:
	Artur Araujo <artur.stat@gmail.com>

Description:
	Computes the mean and variance vectors respectively
		labeled MX and SX.

Parameters:
	T1[in]			pointer to T1 first element.
	SW[in]			pointer to a weights stype structure.
	T2[in]			pointer to T2 first element.
	E[in]			pointer to E first element.
	index[in]		pointer to index first element.
	len[in]			pointer to length of index,
					which must be lower or equal
					than the length of vectors
					T1, SW->ptr, T2, E, MX and SX.
	H[in]			pointer to H first element.
	kfunc[in]		pointer to kernel density function.
	MX[out]			pointer to MX vector.
	SX[out]			pointer to SX vector.
	t[in]			pointer to thread number.
	WORK[out]		pointer to array of transLSW structures.

Return value:
	This function doesn't return a value.

Remarks:
	Vector index must indicate the permutation of vector T2
		sorted by ascending order.
	Vectors T1, SW->ptr, T2, E, MX and SX must have the same length.
*/

static void LSmeasuresI(
	CdoubleCP T1,
	CstypeCP SW,
	CdoubleCP T2,
	CintCP E,
	CintCP index,
	CintCP len,
	Cdouble H[2],
	Kfunc kfunc,
	doubleCP MX,
	doubleCP SX,
	CintCP t,
	transLSW *WORK)
{
	#ifdef _OPENMP
	#pragma omp parallel if( !omp_in_parallel() ) num_threads(global_num_threads)
	#endif
	{
		register int i, j;
		int tid = *t;
		double aux;
		#ifdef _OPENMP
		if (omp_get_num_threads() != 1) tid = omp_get_thread_num(); // use the correct thread number
		#endif
		double *K = WORK[tid].K;
		if (H[1] == H[0]) {
			#ifdef _OPENMP
			#pragma omp for
			#endif
			for (i = 0; i < *len; i++) {
				kfunc(T1, SW, index, &T1[index[i]], &H[0], K); // compute weights
				wkmsurv(len, T2, E, K, index, len, K); // compute conditional survival probabilities vector
				aux = 1-K[index[0]];
				MX[index[i]] = aux*T2[index[0]]; // initialize mean
				SX[index[i]] = aux*R_pow_di(T2[index[0]], 2); // initialize variance
				for (j = 1; j < *len; j++) {
					aux = K[index[j-1]]-K[index[j]]; // compute survival probability jump
					MX[index[i]] += aux*T2[index[j]]; // compute mean
					SX[index[i]] += aux*R_pow_di(T2[index[j]], 2); // compute variance
				}
				aux = 1;
				if (K[index[*len-1]] != 1) aux /= 1-K[index[*len-1]];
				MX[index[i]] *= aux; // normalize mean
				SX[index[i]] *= aux; // normalize variance
				SX[index[i]] -= R_pow_di(MX[index[i]], 2); // compute variance
				if (SX[index[i]] < 1e-10) SX[index[i]] = 1e-10; // variance can't be negative
				SX[index[i]] = sqrt(SX[index[i]]); // compute standard deviation
			}
		} else {
			#ifdef _OPENMP
			#pragma omp for
			#endif
			for (i = 0; i < *len; i++) {
				kfunc(T1, SW, index, &T1[index[i]], &H[0], K); // compute weights
				wkmsurv(len, T2, E, K, index, len, K); // compute conditional survival probabilities vector
				MX[index[i]] = (1-K[index[0]])*T2[index[0]]; // initialize mean
				for (j = 1; j < *len; j++) {
					MX[index[i]] += (K[index[j-1]]-K[index[j]])*T2[index[j]]; // compute mean
				}
				if (K[index[*len-1]] != 1) MX[index[i]] /= 1-K[index[*len-1]]; // normalize mean
			}
			#ifdef _OPENMP
			#pragma omp for
			#endif
			for (i = 0; i < *len; i++) {
				kfunc(T1, SW, index, &T1[index[i]], &H[1], K); // compute weights
				wkmsurv(len, T2, E, K, index, len, K); // compute conditional survival probabilities vector
				SX[index[i]] = (1-K[index[0]])*R_pow_di(T2[index[0]], 2); // initialize variance
				for (j = 1; j < *len; j++) {
					SX[index[i]] += (K[index[j-1]]-K[index[j]])*R_pow_di(T2[index[j]], 2); // compute variance
				}
				if (K[index[*len-1]] != 1) SX[index[i]] /= 1-K[index[*len-1]]; // normalize variance
				SX[index[i]] -= R_pow_di(MX[index[i]], 2); // compute variance
				if (SX[index[i]] < 1e-10) SX[index[i]] = 1e-10; // variance can't be negative
				SX[index[i]] = sqrt(SX[index[i]]); // compute standard deviation
			}
		}
	}
	return;
} // LSmeasuresI

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
	len[in]			pointer to length of T1, E1, T2 and E.
	T1[in]			pointer to T1 first element.
	E1[in]			pointer to E1 first element.
	T2[in]			pointer to T2 first element.
	E[in]			pointer to E first element.
	index0[in]		pointer to index0 first element.
	index1[inout]		pointer to index1 first element.
	nt[in]			pointer to length of UT and number of rows of P.
	UT[in]			pointer to unique times vector.
	nb[in]			pointer to number of rows of P.
	P[out]			pointer to a (nb)x(nt)x4 probability array.
	b[in]			pointer to row index.
	kfunc[in]		pointer to kernel density function.
	H[in]			pointer to H first element.
	t[in]			pointer to thread number.
	WORK[out]		pointer to array of transLSW structures.

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
	int E1[*len],
	Cdouble T2[*len],
	Cint E[*len],
	Cint index0[*len],
	int index1[*len],
	CintCP nt,
	Cdouble UT[*nt],
	CintCP nb,
	double P[*nb*(*nt)*4],
	CintCP b,
	Kfunc kfunc,
	Cdouble H[2],
	CintCP t,
	transLSW *WORK)
{
	register int i;
	int j, k, x, e0, e1;
	double aux[2], p;
	double *W = WORK[*t].a;
	stype SW; // declare stype structure
	SW.type = INT_PTR; // type is an int pointer
	SW.ptr.integer = E1; // hold E1 pointer in ptr union
	SW.length = *len; // hold length of array
	j = 0;
	getIndexI(T1, index0, &UT[0], len, &j, &e0); // determine first index
	for (aux[0] = 1, p = 1; j < e0; j++) { // loop through the sample until last index is reached
		W[index0[j]] = (double)E1[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-W[index0[j]]; // factor needed for the computation
		W[index0[j]] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
		p -= W[index0[j]]; // compute probability
	}
	getIndexI(T1, index0, &UT[*nt-1], len, &j, &e1); // determine last index
	for (i = 0; j < e1; j++) { // loop through the sample until last index is reached
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
	for (; j < *len; j++) { // loop through the sample until last index is reached
		W[index0[j]] = (double)E1[index0[j]]/(*len-j); // compute needed factor
		aux[1] = 1-W[index0[j]]; // factor needed for the computation
		W[index0[j]] *= aux[0]; // compute and save weight
		aux[0] *= aux[1]; // compute and save factor needed for next iteration
	}
	double *MX = WORK[*t].b;
	double *SX = WORK[*t].c;
	LSmeasuresI(T1, &SW, T2, E, index1, len, H, kfunc, MX, SX, t, WORK); // compute location-scale measures
	double *EX = WORK[*t].E0B;
	double *SV = WORK[*t].E1B;
	for (i = 0; i < *len; i++) EX[index1[i]] = (T2[index1[i]]-MX[index1[i]])/SX[index1[i]];
	order_d(EX, index1, *len, FALSE, FALSE, SV); // get permuation
	kmsurv(len, EX, E, index1, len, SV); // compute survival probabilities
	x = e0;
	#ifdef _OPENMP
	#pragma omp parallel for if(*b < 1) num_threads(global_num_threads) private(i, j, k, aux, e1) firstprivate(x) ordered
	#endif
	for (i = 0; i < *nt; i++) {
		#ifdef _OPENMP
		#pragma omp ordered
		#endif
		{
			for (j = 0, P[*b+*nb*(i+*nt*3)] = 0; j < e0; j++) {
				k = 0;
				aux[0] = (UT[i]-T1[index0[j]]-MX[index0[j]])/SX[index0[j]];
				aux[1] = 1;
				getOrdinateI(EX, SV, index1, len, &k, &aux[0], &aux[1]); // get survival probability
				P[*b+*nb*(i+*nt*3)] += aux[1]*W[index0[j]];
			}
			getIndexI(T1, index0, &UT[i], len, &x, &e1); // determine last index
			for (P[*b+*nb*(i+*nt)] = 0; j < e1; j++) {
				k = 0;
				aux[0] = (UT[i]-T1[index0[j]]-MX[index0[j]])/SX[index0[j]];
				aux[1] = 1;
				getOrdinateI(EX, SV, index1, len, &k, &aux[0], &aux[1]); // get survival probability
				P[*b+*nb*(i+*nt)] += aux[1]*W[index0[j]];
			}
			x = e1;
		}
	}
	#ifdef _OPENMP
	#pragma omp parallel for if(*b < 1) num_threads(global_num_threads) private(i) ordered
	#endif
	for (i = *nt-1; i >= 0; i--) {
		#ifdef _OPENMP
		#pragma omp ordered
		#endif
		{
			P[*b+*nb*(i+*nt)] /= P[*b]; // compute and save p12(s,t)
			P[*b+*nb*i] /= P[*b]; // compute and save p11(s,t)
			if (P[*b+*nb*i] < 0) P[*b+*nb*i] = 0;
			P[*b+*nb*(i+*nt*2)] = 1-P[*b+*nb*i]-P[*b+*nb*(i+*nt)]; // compute and save p13(s,t)
			if (P[*b+*nb*(i+*nt*2)] < 0) {
				P[*b+*nb*(i+*nt)] = 1-P[*b+*nb*i]; // compute and save p12(s,t)
				P[*b+*nb*(i+*nt*2)] = 0; // save p13(s,t)
			}
			P[*b+*nb*(i+*nt*3)] /= P[*b+*nb*(*nt*3)]; // compute and save p22(s,t)
		}
	}
	return;
} // transLSI

/*
Author:
	Artur Araujo <artur.stat@gmail.com>

Description:
	Computes a transition probability array
		based on the Location-Scale estimator.

Parameters:
	object			an object of class 'LS'.
	UT			unique times vector.
	h			a vector of bandwidth values with length two.
	nh			number of bandwidth values
					to test by cross-validation.
	ncv			number of cross-validation samples.
	window			a string indicating the desired window or kernel.
	nboot			number of bootstrap samples.
	bootcv			if TRUE cross-validation is done for each bootstrap sample.
	cvfull			if TRUE cross-validation is done for both location and
					scale functions.

Return value:
	Returns a list where the first element is a
		(nboot)x(nt)x4 array of transition probabilities,
		and the second element is a vector of bandwidth
		values used to compute the transition probability
		estimates.
*/

SEXP TransPROBLS(
	SEXP object,
	SEXP UT,
	SEXP h,
	SEXP nh,
	SEXP ncv,
	SEXP window,
	SEXP nboot,
	SEXP bootcv,
	SEXP cvfull)
{
	SEXP data, T1, E1, S, E;
	data = VECTOR_ELT(object, 0);
	T1 = VECTOR_ELT(data, 0);
	E1 = VECTOR_ELT(data, 1);
	S = VECTOR_ELT(data, 2);
	E = VECTOR_ELT(data, 3);
	int len = GET_LENGTH(T1), nt = GET_LENGTH(UT), b, t, nth = 1;
	double *T2 = (double*)malloc( len*sizeof(double) ); // allocate memory block
	if (T2 == NULL) error("TransPROBLS: No more memory\n");
	for (b = 0; b < len; b++) T2[b] = REAL(S)[b]-REAL(T1)[b];
	Kfunc kfunc = kchar2ptr(window); // declare and get pointer to function
	SEXP P, H, list;
	PROTECT( P = alloc3DArray(REALSXP, *INTEGER(nboot), nt, 4) );
	PROTECT( H = NEW_NUMERIC(2) );
	PROTECT( list = NEW_LIST(2) );
	transLSW *WORK = (transLSW*)malloc( global_num_threads*sizeof(transLSW) ); // allocate memory block
	if (WORK == NULL) error("TransPROBLS: No more memory\n");
	for (t = 0; t < global_num_threads; t++) { // allocate per thread memory
		if ( ( WORK[t].sample0 = (int*)malloc( len*sizeof(int) ) ) == NULL ) error("TransPROBLS: No more memory\n");
		if ( ( WORK[t].sample1 = (int*)malloc( len*sizeof(int) ) ) == NULL ) error("TransPROBLS: No more memory\n");
		if ( ( WORK[t].unique0 = (int*)malloc( len*sizeof(int) ) ) == NULL ) error("TransPROBLS: No more memory\n");
		if ( ( WORK[t].K = (double*)malloc( len*sizeof(double) ) ) == NULL ) error("TransPROBLS: No more memory\n");
		if ( ( WORK[t].MX = (double*)malloc( len*sizeof(double) ) ) == NULL ) error("TransPROBLS: No more memory\n");
		if ( ( WORK[t].a = (double*)malloc( len*sizeof(double) ) ) == NULL ) error("TransPROBLS: No more memory\n");
		if ( ( WORK[t].b = (double*)malloc( len*sizeof(double) ) ) == NULL ) error("TransPROBLS: No more memory\n");
		if ( ( WORK[t].c = (double*)malloc( len*sizeof(double) ) ) == NULL ) error("TransPROBLS: No more memory\n");
		if ( ( WORK[t].E0B = (double*)malloc( len*sizeof(double) ) ) == NULL ) error("TransPROBLS: No more memory\n");
		if ( ( WORK[t].E1B = (double*)malloc( len*sizeof(double) ) ) == NULL ) error("TransPROBLS: No more memory\n");
	}
	if (*INTEGER(nboot) > 1) nth = global_num_threads;
	int **index0 = (int**)malloc( nth*sizeof(int*) ); // allocate memory block
	if (index0 == NULL) error("TransPROBLS: No more memory\n");
	int **index1 = (int**)malloc( nth*sizeof(int*) ); // allocate memory block
	if (index1 == NULL) error("TransPROBLS: No more memory\n");
	for (t = 0; t < nth; t++) { // allocate per thread memory
		if ( ( index0[t] = (int*)malloc( len*sizeof(int) ) ) == NULL ) error("TransPROBLS: No more memory\n");
		if ( ( index1[t] = (int*)malloc( len*sizeof(int) ) ) == NULL ) error("TransPROBLS: No more memory\n");
	}
	double HC[2];
	b = 0; // b = len, put it back to 0 or a crash might occur
	t = 0; // t = nth, put it back to 0 or a crash might occur
	indx_ii(&len, index0[0], index1[0]); // initialize indexes
	crossValid(REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), T2, index0[0], &len, REAL(h), INTEGER(nh), INTEGER(ncv), LOGICAL(cvfull), kfunc, HC, &t, WORK); // compute bandwidth
	order_di(REAL(T1), INTEGER(E1), index0[0], len, FALSE, FALSE, TRUE, WORK[0].K, WORK[0].unique0); // get permuation
	order_d(T2, index1[0], len, FALSE, FALSE, WORK[0].K); // get permuation
	transLSI(&len, REAL(T1), INTEGER(E1), T2, INTEGER(E), index0[0], index1[0], &nt, REAL(UT), INTEGER(nboot), REAL(P), &b, kfunc, HC, &t, WORK); // compute transition probabilities
	REAL(H)[0] = HC[0];
	REAL(H)[1] = HC[1];
	if (*INTEGER(nboot) > 1) {
		#ifdef _OPENMP
		#pragma omp parallel num_threads(global_num_threads) private(b, t) firstprivate(HC)
		#endif
		{
			#ifdef _OPENMP
			t = omp_get_thread_num();
			#else
			t = 0;
			#endif
			#ifdef _OPENMP
			#pragma omp for
			#endif
			for (b = 1; b < *INTEGER(nboot); b++) {
				boot_ii(RngArray[t], &len, index0[t], index1[t]); // bootstrap indexes
				if ( *LOGICAL(bootcv) ) crossValid(REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), T2, index0[t], &len, REAL(h), INTEGER(nh), INTEGER(ncv), LOGICAL(cvfull), kfunc, HC, &t, WORK); // compute bandwidth
				order_di(REAL(T1), INTEGER(E1), index0[t], len, FALSE, FALSE, TRUE, WORK[t].K, WORK[t].unique0); // get permuation
				order_d(T2, index1[t], len, FALSE, FALSE, WORK[t].K); // get permuation
				transLSI(&len, REAL(T1), INTEGER(E1), T2, INTEGER(E), index0[t], index1[t], &nt, REAL(UT), INTEGER(nboot), REAL(P), &b, kfunc, HC, &t, WORK); // compute transition probabilities
			}
		}
	}
	for (t = nth-1; t >= 0; t--) {
		free(index0[t]); // free memory block
		free(index1[t]); // free memory block
	}
	free(index0); // free memory block
	free(index1); // free memory block
	for (t = global_num_threads-1; t >= 0; t--) {
		free(WORK[t].sample0); // free memory block
		free(WORK[t].sample1); // free memory block
		free(WORK[t].unique0); // free memory block
		free(WORK[t].K); // free memory block
		free(WORK[t].MX); // free memory block
		free(WORK[t].a); // free memory block
		free(WORK[t].b); // free memory block
		free(WORK[t].c); // free memory block
		free(WORK[t].E0B); // free memory block
		free(WORK[t].E1B); // free memory block
	}
	free(WORK); // free memory block
	free(T2); // free memory block
	SET_ELEMENT(list, 0, P);
	SET_ELEMENT(list, 1, H);
	UNPROTECT(3);
	return list;
} // TransPROBLS
