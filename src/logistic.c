
#include <stdlib.h>
#include <R_ext/Lapack.h>
#include <Rmath.h>
#include "defines.h"

/*
Author:
	Artur Agostinho Araujo <artur.stat@gmail.com>

Description:
	Computes the predicted values.

Parameters:
	len[in]			pointer to length of subset vector.
	subset[in]		pointer to subset vector.
	n[in]			pointer to length of arrays X and B.
	X[in]			pointer to array of pointers.
	B[in]			pointer to vector of parameters.
	P[out]			pointer to predicted values vector.

Return value:
	This function doesn't return a value.

Remarks:
	Vectors X and B must have the same length.
	X must be an array of pointers to the predictor
		variable values.
*/

static void predict(
	CintCP len,
	Cint subset[*len],
	CintCP n,
	doubleCP X[*n],
	Cdouble B[*n],
	doubleCP P)
{
	register int i, j;
	for (i = 0; i < *len; i++) {
		for (P[subset[i]] = 0, j = 0; j < *n; j++) {
			P[subset[i]] += X[j][subset[i]]*B[j]; // compute the linear predictor
		}
		P[subset[i]] = exp(P[subset[i]]);
		P[subset[i]] /= P[subset[i]]+1; // compute the inverse function of the logit
	}
	return;
} // predict

/*
Author:
	Artur Agostinho Araujo <artur.stat@gmail.com>

Description:
	Computes the deviance of the fitted model.

Parameters:
	len[in]			pointer to length of subset vector.
	subset[in]		pointer to subset vector.
	Y[in]			pointer to response values vector.
	P[in]			pointer to predicted values vector.
	d[out]			pointer to deviance value.

Return value:
	This function doesn't return a value.
*/

static void deviance(
	CintCP len,
	Cint subset[*len],
	CintCP Y,
	CdoubleCP P,
	doubleCP d)
{
	register int i;
	for (*d = 0, i = 0; i < *len; i++) {
		if (Y[subset[i]]) *d -= log(P[subset[i]]);
		else *d -= log(1-P[subset[i]]);
	}
	*d *= 2;
	return;
} // deviance

/*
Author:
	Artur Agostinho Araujo <artur.stat@gmail.com>

Description:
	Computes the predicted values from the logistic regression,
	where Y[i]~Bernoulli(P[i]) and logit(P[i])=B[0]+B[1]*X1[i]+B[2]*X2[i]+...

Parameters:
	len[in]			pointer to length of subset vector.
	subset[in]		pointer to subset vector.
	Y[in]			pointer to response values vector.
	P[out]			pointer to predicted values vector.
	n[in]			pointer to length of array X.
	X[in]			pointer to array of pointers.
	maxit[in]		pointer to maximum number of iterations.
	epsilon[in]		pointer to convergence parameter.
	conv[out]		pointer to conv value.

Return value:
	This function doesn't return a value.

Remarks:
	X must be an array of pointers to the predictor
		variable vectors.
	If *conv == 1 the algorithm converged, if *conv == 0 it didn't.
*/

void predict_logit(
	CintCP len,
	Cint subset[*len],
	CintCP Y,
	doubleCP P,
	CintCP n,
	doubleCP X[*n],
	CintCP maxit,
	CdoubleCP epsilon,
	intCP conv)
{
	register int it, i, j, k;
	int info, lwork = (*n)*(*n);
	double dev[2], aux[2];
	double *B = (double*)malloc( *n*sizeof(double) ); // allocate memory block
	double *U = (double*)malloc( *n*sizeof(double) ); // allocate memory block
	double *F = (double*)malloc( lwork*sizeof(double) ); // allocate memory block
	int *IPIV = (int*)malloc( *n*sizeof(int) ); // allocate memory block
	double *WORK = (double*)malloc( lwork*sizeof(double) ); // allocate memory block
	for (i = 0; i < *n; i++) B[i] = 0; // initialize vector of parameters
	predict(len, subset, n, X, B, P); // initialize predicted values
	deviance(len, subset, Y, P, &dev[0]); // initialize deviance
	for (*conv = 0, it = 0; it < *maxit; it++) {
		for (i = 0; i < *n; i++) {
			U[i] = 0; // initialize the score vector
			for (j = 0; j < *n; j++) F[i+*n*j] = 0; // initialize information matrix
		}
		for (k = 0; k < *len; k++) {
			aux[0] = Y[subset[k]]-P[subset[k]];
			aux[1] = P[subset[k]]*(1-P[subset[k]]); // compute variance
			for (i = 0; i < *n; i++) {
				U[i] += X[i][subset[k]]*aux[0]; // compute the score vector
				for (j = i; j < *n; j++) {
					F[i+*n*j] += X[i][subset[k]]*X[j][subset[k]]*aux[1]; // fill upper triangular part of the information matrix
				}
			}
		}
		for (i = 0; i < *n-1; i++) {
			for (j = i+1; j < *n; j++) {
				F[j+*n*i] = F[i+*n*j]; // copy upper part of the information matrix to it's lower part
			}
		}
		F77_CALL(dgetrf)(n, n, F, n, IPIV, &info); // compute LU factorization of the information matrix
		F77_CALL(dgetri)(n, F, n, IPIV, WORK, &lwork, &info); // compute the inverse of the information matrix
		for (i = 0; i < *n; i++) {
			for (j = 0; j < *n; j++) {
				B[i] += F[i+*n*j]*U[j]; // compute vector of parameters
			}
		}
		predict(len, subset, n, X, B, P); // compute predicted values
		deviance(len, subset, Y, P, &dev[1]); // compute deviance
		if (isnan(dev[1]) || info) { // if isnan(dev[1]) the algorithm isn't converging
			for (i = 0; i < *len; i++) P[subset[i]] = Y[subset[i]]; // in that case make the predicted values equal to the response values
			break;
		}
		if (fabs(dev[1]-dev[0])/( 0.1 + fabs(dev[1]) ) < *epsilon) { // check for convergence
			*conv = 1;
			break;
		} else dev[0] = dev[1];
	}
	free(B); // free memory block
	free(U); // free memory block
	free(F); // free memory block
	free(IPIV); // free memory block
	free(WORK); // free memory block
	return;
} // predict_logit
