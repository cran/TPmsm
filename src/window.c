
#include <stdlib.h>
#include <string.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "defines.h"
#include "sort.h"

#ifndef M_70_81
#define M_70_81	0.864197530864197482891597701382	/* 70/81 */
#endif

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes weights based on the gaussian kernel.

Parameters:
	X[in]			pointer to covariate vector.
	n[in]			pointer to lenght of covariate vector.
	x[in]			pointer to covariate value to compute the weights at.
	h[in]			pointer to bandwith value.
	W[out]			pointer to weights vector.

Return value:
	This function doesn't return a value.
*/

static void knormal(
	CdoubleCP X,
	CintCP n,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP W)
{
	register int i;
	for (i = 0; i < *n; i++) {
		W[i] = (X[i] - *x) / *h;
		W[i] = M_1_SQRT_2PI*exp( -0.5*R_pow_di(W[i], 2) ); // 1/sqrt(2pi) = M_1_SQRT_2PI, 1/2 = 0.5
	}
	return;
} // knormal

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes weights based on the epanechnikov kernel.

Parameters:
	X[in]			pointer to covariate vector.
	n[in]			pointer to lenght of covariate vector.
	x[in]			pointer to covariate value to compute the weights at.
	h[in]			pointer to bandwith value.
	W[out]			pointer to weights vector.

Return value:
	This function doesn't return a value.
*/

static void kepanech(
	CdoubleCP X,
	CintCP n,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP W)
{
	register int i;
	for (i = 0; i < *n; i++) {
		W[i] = (X[i] - *x) / *h;
		W[i] = 0.75*( 1-R_pow_di(W[i], 2) )*(fabs(W[i]) <= 1); // 3/4 = 0.75
	}
	return;
} // kepanech

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes weights based on the tricube kernel.

Parameters:
	X[in]			pointer to covariate vector.
	n[in]			pointer to lenght of covariate vector.
	x[in]			pointer to covariate value to compute the weights at.
	h[in]			pointer to bandwith value.
	W[out]			pointer to weights vector.

Return value:
	This function doesn't return a value.
*/

static void ktricube(
	CdoubleCP X,
	CintCP n,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP W)
{
	register int i;
	for (i = 0; i < *n; i++) {
		W[i] = (X[i] - *x) / *h;
		W[i] = fabs(W[i]);
		W[i] = M_70_81*R_pow_di(1-R_pow_di(W[i], 3), 3)*(W[i] <= 1); // 70/81 = M_70_81
	}
	return;
} // ktricube

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes weights based on the boxcar kernel.

Parameters:
	X[in]			pointer to covariate vector.
	n[in]			pointer to lenght of covariate vector.
	x[in]			pointer to covariate value to compute the weights at.
	h[in]			pointer to bandwith value.
	W[out]			pointer to weights vector.

Return value:
	This function doesn't return a value.
*/

static void kbox(
	CdoubleCP X,
	CintCP n,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP W)
{
	register int i;
	for (i = 0; i < *n; i++) {
		W[i] = (X[i] - *x) / *h;
		W[i] = 0.5*(fabs(W[i]) <= 1); // 1/2 = 0.5
	}
	return;
} // kbox

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes weights based on the triangular kernel.

Parameters:
	X[in]			pointer to covariate vector.
	n[in]			pointer to lenght of covariate vector.
	x[in]			pointer to covariate value to compute the weights at.
	h[in]			pointer to bandwith value.
	W[out]			pointer to weights vector.

Return value:
	This function doesn't return a value.
*/

static void ktriangular(
	CdoubleCP X,
	CintCP n,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP W)
{
	register int i;
	for (i = 0; i < *n; i++) {
		W[i] = (X[i] - *x) / *h;
		W[i] = fabs(W[i]);
		W[i] = (1-W[i])*(W[i] <= 1);
	}
	return;
} // ktriangular

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes weights based on the biweight kernel.

Parameters:
	X[in]			pointer to covariate vector.
	n[in]			pointer to lenght of covariate vector.
	x[in]			pointer to covariate value to compute the weights at.
	h[in]			pointer to bandwith value.
	W[out]			pointer to weights vector.

Return value:
	This function doesn't return a value.
*/

static void kbiweight(
	CdoubleCP X,
	CintCP n,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP W)
{
	register int i;
	for (i = 0; i < *n; i++) {
		W[i] = (X[i] - *x) / *h;
		W[i] = 0.9375*R_pow_di(1-R_pow_di(W[i], 2), 2)*(fabs(W[i]) <= 1); // 15/16 = 0.9375
	}
	return;
} // kbiweight

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes weights based on the triweight kernel.

Parameters:
	X[in]			pointer to covariate vector.
	n[in]			pointer to lenght of covariate vector.
	x[in]			pointer to covariate value to compute the weights at.
	h[in]			pointer to bandwith value.
	W[out]			pointer to weights vector.

Return value:
	This function doesn't return a value.
*/

static void ktriweight(
	CdoubleCP X,
	CintCP n,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP W)
{
	register int i;
	for (i = 0; i < *n; i++) {
		W[i] = (X[i] - *x) / *h;
		W[i] = 1.09375*R_pow_di(1-R_pow_di(W[i], 2), 3)*(fabs(W[i]) <= 1); // 35/32 = 1.09375
	}
	return;
} // ktriweight

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes weights based on the cosine kernel.

Parameters:
	X[in]			pointer to covariate vector.
	n[in]			pointer to lenght of covariate vector.
	x[in]			pointer to covariate value to compute the weights at.
	h[in]			pointer to bandwith value.
	W[out]			pointer to weights vector.

Return value:
	This function doesn't return a value.
*/

static void kcosine(
	CdoubleCP X,
	CintCP n,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP W)
{
	register int i;
	for (i = 0; i < *n; i++) {
		W[i] = (X[i] - *x) / *h;
		W[i] = M_PI_4*cos(M_PI_2*W[i])*(fabs(W[i]) <= 1); // pi/4 = M_PI_4, pi/2 = M_PI_2
	}
	return;
} // kcosine

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes nearest neighbour weights based on an asymmetric window.

Parameters:
	X[in]			pointer to covariate vector.
	n[in]			pointer to lenght of covariate vector.
	x[in]			pointer to covariate value to compute the weights at.
	h[in]			pointer to bandwith value.
	W[out]			pointer to weights vector.

Return value:
	This function doesn't return a value.
*/

static void wasymmetric(
	CdoubleCP X,
	CintCP n,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP W)
{
	register int i;
	int *index = (int*)malloc( *n*sizeof(int) ), tr;
	double *WORK = (double*)malloc( *n*sizeof(double) ); // allocate memory block
	double lambda;
	for (i = 0; i < *n; i++) {
		W[i] = X[i] - *x; // initialize weights vector
		index[i] = i; // initialize index vector
	}
	order_d(W, index, *n, FALSE, FALSE, WORK); // get permutation
	free(WORK); // free memory block
	tr = *n * *h;
	for (i = *n-1; i > 0; i--) {
		if (W[index[i]] < 0) break; // compute index
	}
	if (i+tr+1 > *n-1) lambda = W[index[*n-1]];
	else lambda = W[index[i+tr+1]];
	for (i = 0; i < *n; i++) {
		W[index[i]] = W[index[i]] >= 0 && W[index[i]] <= lambda; // compute weights
	}
	free(index); // free memory block
	return;
} // wasymmetric

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes nearest neighbour weights based on a symmetric window.

Parameters:
	X[in]			pointer to covariate vector.
	n[in]			pointer to lenght of covariate vector.
	x[in]			pointer to covariate value to compute the weights at.
	h[in]			pointer to bandwith value.
	W[out]			pointer to weights vector.

Return value:
	This function doesn't return a value.
*/

static void wsymmetric(
	CdoubleCP X,
	CintCP n,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP W)
{
	register int i;
	int *index = (int*)malloc( *n*sizeof(int) ), tr;
	double *WORK = (double*)malloc( *n*sizeof(double) ); // allocate memory block
	double lambda[2];
	for (i = 0; i < *n; i++) {
		W[i] = X[i] - *x; // initialize weights vector
		index[i] = i; // initialize index vector
	}
	order_d(W, index, *n, FALSE, FALSE, WORK); // get permutation
	free(WORK); // free memory block
	tr = *n * *h / 2;
	for (i = *n-1; i > 0; i--) {
		if (W[index[i]] <= 0) break; // compute lower index
	}
	if (i-tr < 0) lambda[0] = -fabs(W[index[0]]);
	else lambda[0] = -fabs(W[index[i-tr]]);
	for (; i > 0; i--) {
		if (W[index[i]] < 0) break; // compute upper index
	}
	if (i+tr+1 > *n-1) lambda[1] = W[index[*n-1]];
	else lambda[1] = W[index[i+tr+1]];
	for (i = 0; i < *n; i++) {
		W[index[i]] = W[index[i]] >= lambda[0] && W[index[i]] <= lambda[1]; // compute weights
	}
	free(index); // free memory block
	return;
} // wsymmetric

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Returns a pointer to a window or kernel function based on the inputed string.

Parameters:
	window			a string indicating the desired window or kernel method.

Return value:
	Returns a pointer to a window or kernel function.
*/

void (*kchar2ptr(SEXP window))(CdoubleCP, CintCP, CdoubleCP, CdoubleCP, doubleCP) {
	CcharCP window1 = CHAR( STRING_ELT(window, 0) );
	if (strcmp(window1, "epanech") == 0) return kepanech; // if window is epanechnikov
	else if (strcmp(window1, "tricube") == 0) return ktricube; // if window is tricube
	else if (strcmp(window1, "box") == 0) return kbox; // if window is boxcar
	else if (strcmp(window1, "triangular") == 0) return ktriangular; // if window is triangular
	else if (strcmp(window1, "biweight") == 0) return kbiweight; // if window is biweight
	else if (strcmp(window1, "triweight") == 0) return ktriweight; // if window is triweight
	else if (strcmp(window1, "cosine") == 0) return kcosine; // if window is cosine
	else if (strcmp(window1, "asymmetric") == 0) return wasymmetric; // if window is asymmetric
	else if (strcmp(window1, "symmetric") == 0) return wsymmetric; // if window is symmetric
	else return knormal; // defaults to gaussian
} // kchar2ptr

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes the Nadaraya-Watson weights.

Parameters:
	X[in]			pointer to covariate vector.
	n[in]			pointer to lenght of covariate vector.
	x[in]			pointer to covariate value to compute the weights at.
	h[in]			pointer to bandwith value.
	W[out]			pointer to weights vector.
	kfunc[in] 		pointer to kernel density function.

Return value:
	This function doesn't return a value.
*/

void NWWeights(
	CdoubleCP X,
	CintCP n,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP W,
	void (*kfunc)(CdoubleCP, CintCP, CdoubleCP, CdoubleCP, doubleCP) )
{
	register int i;
	double sum;
	kfunc(X, n, x, h, W); // compute kernel density
	for (i = 0, sum = 0; i < *n; i++) sum += W[i]; // compute sum
	for (i = 0; i < *n; i++) W[i] /= sum; // compute Nadaraya-Watson weight
	return;
} // NWWeights

/*
Author:
	Artur Agostinho Araújo <b5498@math.uminho.pt>

Description:
	Computes local linear weights based on a kernel density.

Parameters:
	X[in]			pointer to covariate vector.
	n[in]			pointer to lenght of covariate vector.
	x[in]			pointer to covariate value to compute the weights at.
	h[in]			pointer to bandwith value.
	W[out]			pointer to weights vector.
	kfunc[in] 		pointer to kernel density function.

Return value:
	This function doesn't return a value.
*/

void LLWeights(
	CdoubleCP X,
	CintCP n,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP W,
	void (*kfunc)(CdoubleCP, CintCP, CdoubleCP, CdoubleCP, doubleCP) )
{
	register int i;
	double aux[2], sum[2] = {0, 0};
	kfunc(X, n, x, h, W); // compute kernel density
	for (i = 0; i < *n; i++) {
		aux[0] = X[i]-*x;
		aux[1] = aux[0]*W[i];
		sum[0] += aux[1];
		sum[1] += aux[0]*aux[1];
	}
	for (i = 0; i < *n; i++) W[i] *= sum[1]-(X[i]-*x)*sum[0];
	for (i = 0, sum[0] = 0; i < *n; i++) sum[0] += W[i]; // compute sum
	for (i = 0; i < *n; i++) W[i] /= sum[0]; // compute local linear weight
	return;
} // LLWeights
