
#include <string.h>
#include <Rdefines.h>
#include <Rmath.h>
#include "defines.h"
#include "wtypefunc.h"

#ifndef M_70_81
#define M_70_81	0.864197530864197482891597701382	/* 70/81 */
#endif

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Initialize kernel density with weights over bandwidth value.

Parameters:
  SW[in]            pointer to a weights stype structure.
  index[i]          pointer to index vector
  h[in]             pointer to bandwidth value.
  K[out]            pointer to kernel density vector.

Return value:
  This function doesn't return a value.

Remarks:
  The length of the vectors is read from SW->length.
*/

static void kweight(
	CstypeCP SW,
	CintCP index,
	CdoubleCP h,
	doubleCP K)
{
	register int i;
	switch (SW->type) {
		case SINT_PTR:
			K[index[0]] = SW->ptr.shortinteger[index[0]] / *h;
			for (i = 1; i < SW->length; i++) {
				if (index[i] == index[i-1]) continue; // write only once to an index
				K[index[i]] = SW->ptr.shortinteger[index[i]] / *h;
			}
			break;
		case INT_PTR:
			K[index[0]] = SW->ptr.integer[index[0]] / *h;
			for (i = 1; i < SW->length; i++) {
				if (index[i] == index[i-1]) continue; // write only once to an index
				K[index[i]] = SW->ptr.integer[index[i]] / *h;
			}
			break;
		case REAL_PTR:
			K[index[0]] = SW->ptr.real[index[0]] / *h;
			for (i = 1; i < SW->length; i++) {
				if (index[i] == index[i-1]) continue; // write only once to an index
				K[index[i]] = SW->ptr.real[index[i]] / *h;
			}
			break;
		default:
			break;
	}
	return;
} // kweight

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Computes weights based on the gaussian kernel.

Parameters:
  X[in]             pointer to covariate vector.
  SW[in]            pointer to a weights stype structure.
  index[i]          pointer to index vector
  x[in]             pointer to covariate value to compute the weights at.
  h[in]             pointer to bandwidth value.
  K[out]            pointer to kernel density vector.

Return value:
  This function doesn't return a value.

Remarks:
  The length of the vectors is read from SW->length.
*/

static void knormal(
	CdoubleCP X,
	CstypeCP SW,
	CintCP index,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP K)
{
	register int i;
	kweight(SW, index, h, K); // initialize kernel density
	K[index[0]] *= M_1_SQRT_2PI*exp( -0.5*R_pow_di( (X[index[0]] - *x) / *h, 2 ) ); // 1/sqrt(2pi) = M_1_SQRT_2PI, 1/2 = 0.5
	for (i = 1; i < SW->length; i++) {
		if (index[i] == index[i-1]) continue; // write only once to an index
		K[index[i]] *= M_1_SQRT_2PI*exp( -0.5*R_pow_di( (X[index[i]] - *x) / *h, 2 ) ); // 1/sqrt(2pi) = M_1_SQRT_2PI, 1/2 = 0.5
	}
	return;
} // knormal

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Computes weights based on the epanechnikov kernel.

Parameters:
  X[in]             pointer to covariate vector.
  SW[in]            pointer to a weights stype structure.
  index[i]          pointer to index vector
  x[in]             pointer to covariate value to compute the weights at.
  h[in]             pointer to bandwidth value.
  K[out]            pointer to kernel density vector.

Return value:
  This function doesn't return a value.

Remarks:
  The length of the vectors is read from SW->length.
*/

static void kepanech(
	CdoubleCP X,
	CstypeCP SW,
	CintCP index,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP K)
{
	register int i;
	double aux;
	kweight(SW, index, h, K); // initialize kernel density
	aux = (X[index[0]] - *x) / *h;
	K[index[0]] *= 0.75*( 1-R_pow_di(aux, 2) )*(fabs(aux) <= 1); // 3/4 = 0.75
	for (i = 1; i < SW->length; i++) {
		if (index[i] == index[i-1]) continue; // write only once to an index
		aux = (X[index[i]] - *x) / *h;
		K[index[i]] *= 0.75*( 1-R_pow_di(aux, 2) )*(fabs(aux) <= 1); // 3/4 = 0.75
	}
	return;
} // kepanech

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Computes weights based on the tricube kernel.

Parameters:
  X[in]             pointer to covariate vector.
  SW[in]            pointer to a weights stype structure.
  index[i]          pointer to index vector
  x[in]             pointer to covariate value to compute the weights at.
  h[in]             pointer to bandwidth value.
  K[out]            pointer to kernel density vector.

Return value:
  This function doesn't return a value.

Remarks:
  The length of the vectors is read from SW->length.
*/

static void ktricube(
	CdoubleCP X,
	CstypeCP SW,
	CintCP index,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP K)
{
	register int i;
	double aux;
	kweight(SW, index, h, K); // initialize kernel density
	aux = fabs( (X[index[0]] - *x) / *h );
	K[index[0]] *= M_70_81*R_pow_di(1-R_pow_di(aux, 3), 3)*(aux <= 1); // 70/81 = M_70_81
	for (i = 1; i < SW->length; i++) {
		if (index[i] == index[i-1]) continue; // write only once to an index
		aux = fabs( (X[index[i]] - *x) / *h );
		K[index[i]] *= M_70_81*R_pow_di(1-R_pow_di(aux, 3), 3)*(aux <= 1); // 70/81 = M_70_81
	}
	return;
} // ktricube

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Computes weights based on the boxcar kernel.

Parameters:
  X[in]             pointer to covariate vector.
  SW[in]            pointer to a weights stype structure.
  index[i]          pointer to index vector
  x[in]             pointer to covariate value to compute the weights at.
  h[in]             pointer to bandwidth value.
  K[out]            pointer to kernel density vector.

Return value:
  This function doesn't return a value.

Remarks:
  The length of the vectors is read from SW->length.
*/

static void kbox(
	CdoubleCP X,
	CstypeCP SW,
	CintCP index,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP K)
{
	register int i;
	kweight(SW, index, h, K); // initialize kernel density
	K[index[0]] *= 0.5*(fabs( (X[index[0]] - *x) / *h ) <= 1); // 1/2 = 0.5
	for (i = 1; i < SW->length; i++) {
		if (index[i] == index[i-1]) continue; // write only once to an index
		K[index[i]] *= 0.5*(fabs( (X[index[i]] - *x) / *h ) <= 1); // 1/2 = 0.5
	}
	return;
} // kbox

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Computes weights based on the triangular kernel.

Parameters:
  X[in]             pointer to covariate vector.
  SW[in]            pointer to a weights stype structure.
  index[i]          pointer to index vector
  x[in]             pointer to covariate value to compute the weights at.
  h[in]             pointer to bandwidth value.
  K[out]            pointer to kernel density vector.

Return value:
  This function doesn't return a value.

Remarks:
  The length of the vectors is read from SW->length.
*/

static void ktriangular(
	CdoubleCP X,
	CstypeCP SW,
	CintCP index,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP K)
{
	register int i;
	double aux;
	kweight(SW, index, h, K); // initialize kernel density
	aux = fabs( (X[index[0]] - *x) / *h );
	K[index[0]] *= (1-aux)*(aux <= 1);
	for (i = 1; i < SW->length; i++) {
		if (index[i] == index[i-1]) continue; // write only once to an index
		aux = fabs( (X[index[i]] - *x) / *h );
		K[index[i]] *= (1-aux)*(aux <= 1);
	}
	return;
} // ktriangular

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Computes weights based on the biweight kernel.

Parameters:
  X[in]             pointer to covariate vector.
  SW[in]            pointer to a weights stype structure.
  index[i]          pointer to index vector
  x[in]             pointer to covariate value to compute the weights at.
  h[in]             pointer to bandwidth value.
  K[out]            pointer to kernel density vector.

Return value:
  This function doesn't return a value.

Remarks:
  The length of the vectors is read from SW->length.
*/

static void kbiweight(
	CdoubleCP X,
	CstypeCP SW,
	CintCP index,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP K)
{
	register int i;
	double aux;
	kweight(SW, index, h, K); // initialize kernel density
	aux = (X[index[0]] - *x) / *h;
	K[index[0]] *= 0.9375*R_pow_di(1-R_pow_di(aux, 2), 2)*(fabs(aux) <= 1); // 15/16 = 0.9375
	for (i = 1; i < SW->length; i++) {
		if (index[i] == index[i-1]) continue; // write only once to an index
		aux = (X[index[i]] - *x) / *h;
		K[index[i]] *= 0.9375*R_pow_di(1-R_pow_di(aux, 2), 2)*(fabs(aux) <= 1); // 15/16 = 0.9375
	}
	return;
} // kbiweight

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Computes weights based on the triweight kernel.

Parameters:
  X[in]             pointer to covariate vector.
  SW[in]            pointer to a weights stype structure.
  index[i]          pointer to index vector
  x[in]             pointer to covariate value to compute the weights at.
  h[in]             pointer to bandwidth value.
  K[out]            pointer to kernel density vector.

Return value:
  This function doesn't return a value.

Remarks:
  The length of the vectors is read from SW->length.
*/

static void ktriweight(
	CdoubleCP X,
	CstypeCP SW,
	CintCP index,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP K)
{
	register int i;
	double aux;
	kweight(SW, index, h, K); // initialize kernel density
	aux = (X[index[0]] - *x) / *h;
	K[index[0]] *= 1.09375*R_pow_di(1-R_pow_di(aux, 2), 3)*(fabs(aux) <= 1); // 35/32 = 1.09375
	for (i = 1; i < SW->length; i++) {
		if (index[i] == index[i-1]) continue; // write only once to an index
		aux = (X[index[i]] - *x) / *h;
		K[index[i]] *= 1.09375*R_pow_di(1-R_pow_di(aux, 2), 3)*(fabs(aux) <= 1); // 35/32 = 1.09375
	}
	return;
} // ktriweight

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Computes weights based on the cosine kernel.

Parameters:
  X[in]             pointer to covariate vector.
  SW[in]            pointer to a weights stype structure.
  index[i]          pointer to index vector
  x[in]             pointer to covariate value to compute the weights at.
  h[in]             pointer to bandwidth value.
  K[out]            pointer to kernel density vector.

Return value:
  This function doesn't return a value.

Remarks:
  The length of the vectors is read from SW->length.
*/

static void kcosine(
	CdoubleCP X,
	CstypeCP SW,
	CintCP index,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP K)
{
	register int i;
	double aux;
	kweight(SW, index, h, K); // initialize kernel density
	aux = (X[index[0]] - *x) / *h;
	K[index[0]] *= M_PI_4*cos(M_PI_2*aux)*(fabs(aux) <= 1); // pi/4 = M_PI_4, pi/2 = M_PI_2
	for (i = 1; i < SW->length; i++) {
		if (index[i] == index[i-1]) continue; // write only once to an index
		aux = (X[index[i]] - *x) / *h;
		K[index[i]] *= M_PI_4*cos(M_PI_2*aux)*(fabs(aux) <= 1); // pi/4 = M_PI_4, pi/2 = M_PI_2
	}
	return;
} // kcosine

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Returns a pointer to a window or kernel function based on the inputted string.

Parameters:
  window            a string indicating the desired window or kernel method.

Return value:
  Returns a pointer to a window or kernel function.
*/

void (*kchar2ptr(SEXP window))(CdoubleCP, CstypeCP, CintCP, CdoubleCP, CdoubleCP, doubleCP) {
	CcharCP window1 = CHAR( STRING_ELT(window, 0) );
	if (strcmp(window1, "epanech") == 0) return kepanech; // if window is epanechnikov
	else if (strcmp(window1, "tricube") == 0) return ktricube; // if window is tricube
	else if (strcmp(window1, "box") == 0) return kbox; // if window is boxcar
	else if (strcmp(window1, "triangular") == 0) return ktriangular; // if window is triangular
	else if (strcmp(window1, "biweight") == 0) return kbiweight; // if window is biweight
	else if (strcmp(window1, "triweight") == 0) return ktriweight; // if window is triweight
	else if (strcmp(window1, "cosine") == 0) return kcosine; // if window is cosine
	else return knormal; // defaults to gaussian
} // kchar2ptr

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Computes the Nadaraya-Watson weights.

Parameters:
  X[in]             pointer to covariate vector.
  SW[in]            pointer to a weights stype structure.
  index[i]          pointer to index vector
  x[in]             pointer to covariate value to compute the weights at.
  h[in]             pointer to bandwidth value.
  K[out]            pointer to kernel density vector.
  kfunc[in]         pointer to kernel density function.

Return value:
  This function doesn't return a value.

Remarks:
  Vector index must indicate the permutation of vector X
    sorted by ascending order.
  Actually index can indicate the permutation of another
    vector sorted by ascending order. As long as that
    vector is indexed to vector X. As the columns in
    a data.frame or matrix are indexed to each other.
  Vectors X, SW->ptr, index and K must have the same length.
  The length of the vectors is read from SW->length.
*/

void NWWeights(
	CdoubleCP X,
	CstypeCP SW,
	CintCP index,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP K,
	Kfunc kfunc)
{
	register int i;
	double sum;
	kfunc(X, SW, index, x, h, K); // compute kernel density
	for (i = 0, sum = 0; i < SW->length; i++) sum += K[index[i]]; // compute sum
	for (K[index[0]] /= sum, i = 1; i < SW->length; i++) {
		if (index[i] == index[i-1]) continue; // write only once to an index
		K[index[i]] /= sum; // compute Nadaraya-Watson weight
	}
	return;
} // NWWeights

/*
Author:
  Artur Araujo <artur.stat@gmail.com>

Description:
  Computes local linear weights based on a kernel density.

Parameters:
  X[in]             pointer to covariate vector.
  SW[in]            pointer to a weights stype structure.
  index[i]          pointer to index vector
  x[in]             pointer to covariate value to compute the weights at.
  h[in]             pointer to bandwidth value.
  K[out]            pointer to kernel density vector.
  kfunc[in]         pointer to kernel density function.

Return value:
  This function doesn't return a value.

Remarks:
  Vector index must indicate the permutation of vector X
    sorted by ascending order.
  Actually index can indicate the permutation of another
    vector sorted by ascending order. As long as that
    vector is indexed to vector X. As the columns in
    a data.frame or matrix are indexed to each other.
  Vectors X, SW->ptr, index and K must have the same length.
  The length of the vectors is read from SW->length.
*/

void LLWeights(
	CdoubleCP X,
	CstypeCP SW,
	CintCP index,
	CdoubleCP x,
	CdoubleCP h,
	doubleCP K,
	Kfunc kfunc)
{
	register int i;
	double aux[2], sum[2] = {0, 0};
	kfunc(X, SW, index, x, h, K); // compute kernel density
	for (i = 0; i < SW->length; i++) {
		aux[0] = X[index[i]] - *x;
		aux[1] = aux[0]*K[index[i]];
		sum[0] += aux[1];
		sum[1] += aux[0]*aux[1];
	}
	K[index[0]] *= sum[1]-(X[index[0]] - *x)*sum[0];
	for (i = 1; i < SW->length; i++) {
		if (index[i] == index[i-1]) continue; // write only once to an index
		K[index[i]] *= sum[1]-(X[index[i]] - *x)*sum[0];
	}
	for (i = 0, sum[0] = 0; i < SW->length; i++) sum[0] += K[index[i]]; // compute sum
	for (K[index[0]] /= sum[0], i = 1; i < SW->length; i++) {
		if (index[i] == index[i-1]) continue; // write only once to an index
		K[index[i]] /= sum[0]; // compute local linear weight
	}
	return;
} // LLWeights
