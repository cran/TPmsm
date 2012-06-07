
#include <errno.h>
#include "defines.h"

/*
 *	Natural Splines
 *	---------------
 *	Here the end-conditions are determined by setting the second
 *	derivative of the spline at the end-points to equal to zero.
 *
 *	There are *n-2 unknowns (y[index0[i]]'' at x[index0[1]], ..., x[index0[*n-2]])
 *  and *n-2 equations to determine them.  Either Choleski or Gaussian
 *	elimination could be used.
 */

static void natural_splineI(
	CdoubleCP x,
	CdoubleCP y,
	CintCP index0,
	CintCP n,
	double b[*n],
	double c[*n],
	double d[*n])
{
	int nm1, i;
	double t;

	if (*n < 2) {
		errno = EDOM;
		return;
    }

	if (*n < 3) {
		t = (y[index0[1]]-y[index0[0]]);
		b[0] = t/(x[index0[1]]-x[index0[0]]);
		b[1] = b[0];
		c[0] = c[1] = d[0] = d[1] = 0.0;
		return;
    }

	nm1 = *n-2;

	/* Set up the tridiagonal system */
	/* b = diagonal, d = offdiagonal, c = right hand side */

	d[0] = x[index0[1]]-x[index0[0]];
	c[1] = (y[index0[1]]-y[index0[0]])/d[0];
	for (i = 1; i < *n-1; i++) {
		d[i] = x[index0[i+1]]-x[index0[i]];
		b[i] = 2.0*(d[i-1]+d[i]);
		c[i+1] = (y[index0[i+1]]-y[index0[i]])/d[i];
		c[i] = c[i+1]-c[i];
	}

	/* Gaussian elimination */

	for (i = 2; i < *n-1; i++) {
		t = d[i-1]/b[i-1];
		b[i] = b[i]-t*d[i-1];
		c[i] = c[i]-t*c[i-1];
	}

	/* Backward substitution */

	c[nm1] = c[nm1]/b[nm1];
	for (i = *n-3; i > 0; i--) c[i] = (c[i]-d[i]*c[i+1])/b[i];

	/* End conditions */

	c[0] = c[*n-1] = 0.0;

	/* Get cubic coefficients */

	b[0] = (y[index0[1]]-y[index0[0]])/d[0]-d[i]*c[1];
	c[0] = 0.0;
	d[0] = c[1]/d[0];
	b[*n-1] = (y[index0[*n-1]]-y[index0[nm1]])/d[nm1]+d[nm1]*c[nm1];
	for (i = 1; i < *n-1; i++) {
		b[i] = (y[index0[i+1]]-y[index0[i]])/d[i]-d[i]*(c[i+1]+2.0*c[i]);
		d[i] = (c[i+1]-c[i])/d[i];
		c[i] = 3.0*c[i];
	}
	c[*n-1] = 0.0;
	d[*n-1] = 0.0;

	return;
} // natural_splineI

/*
 *	Splines a la Forsythe Malcolm and Moler
 *	---------------------------------------
 *	In this case the end-conditions are determined by fitting
 *	cubic polynomials to the first and last 4 points and matching
 *	the third derivitives of the spline at the end-points to the
 *	third derivatives of these cubics at the end-points.
 */

static void fmm_splineI(
	CdoubleCP x,
	CdoubleCP y,
	CintCP index0,
	CintCP n,
	double b[*n],
	double c[*n],
	double d[*n])
{
	int nm1, i;
	double t;

    if (*n < 2) {
		errno = EDOM;
		return;
    }

	if (*n < 3) {
		t = (y[index0[1]]-y[index0[0]]);
		b[0] = t/(x[index0[1]]-x[index0[0]]);
		b[1] = b[0];
		c[0] = c[1] = d[0] = d[1] = 0.0;
		return;
	}

	nm1 = *n-2;

	/* Set up tridiagonal system */
	/* b = diagonal, d = offdiagonal, c = right hand side */

	d[0] = x[index0[1]]-x[index0[0]];
	c[1] = (y[index0[1]]-y[index0[0]])/d[0]; /* = +/- Inf	for x[index0[0]]=x[index0[1]] -- problem? */
	for (i = 1; i < *n-1; i++) {
		d[i] = x[index0[i+1]]-x[index0[i]];
		b[i] = 2.0*(d[i-1]+d[i]);
		c[i+1] = (y[index0[i+1]]-y[index0[i]])/d[i];
		c[i] = c[i+1]-c[i];
	}

	/* End conditions. */
	/* Third derivatives at x[index0[0]] and x[index0[*n-1]] obtained */
	/* from divided differences */

	b[0] = -d[0];
	b[*n-1] = -d[nm1];
	c[0] = c[*n-1] = 0.0;
	if (*n > 3) {
		c[0] = c[2]/(x[index0[3]]-x[index0[1]])-c[1]/(x[index0[2]]-x[index0[0]]);
		c[*n-1] = c[nm1]/(x[index0[*n-1]]-x[index0[*n-3]])-c[*n-3]/(x[index0[nm1]]-x[index0[*n-4]]);
		c[0] = c[0]*d[0]*d[0]/(x[index0[3]]-x[index0[0]]);
		c[*n-1] = -c[*n-1]*d[nm1]*d[nm1]/(x[index0[*n-1]]-x[index0[*n-4]]);
    }

	/* Gaussian elimination */

	for (i = 1; i <= *n-1; i++) {
		t = d[i-1]/b[i-1];
		b[i] = b[i]-t*d[i-1];
		c[i] = c[i]-t*c[i-1];
	}

	/* Backward substitution */

	c[*n-1] = c[*n-1]/b[*n-1];
	for (i = nm1; i >= 0; i--) c[i] = (c[i]-d[i]*c[i+1])/b[i];

	/* c[i] is now the sigma[i-1] of the text */
	/* Compute polynomial coefficients */

    b[*n-1] = (y[index0[*n-1]]-y[index0[*n-2]])/d[*n-2]+d[*n-2]*(c[*n-2]+2.0*c[*n-1]);
    for (i = 0; i <= nm1; i++) {
		b[i] = (y[index0[i+1]]-y[index0[i]])/d[i]-d[i]*(c[i+1]+2.0*c[i]);
		d[i] = (c[i+1]-c[i])/d[i];
		c[i] = 3.0*c[i];
    }
    c[*n-1] = 3.0*c[*n-1];
    d[*n-1] = d[nm1];

	return;
} // fmm_splineI

void spline_coefI(
	CintCP method,
	CdoubleCP x,
	CdoubleCP y,
	CintCP index0,
	CintCP n,
	double b[*n],
	double c[*n],
	double d[*n])
{
	switch (*method) {
		case 0:
			natural_splineI(x, y, index0, n, b, c, d);
			break;
		default:
			fmm_splineI(x, y, index0, n, b, c, d);
	}
	return;
} // spline_coefI

void spline_evalI(
	CintCP method,
	CdoubleCP x,
	CdoubleCP y,
	CintCP index0,
	CintCP n,
	Cdouble b[*n],
	Cdouble c[*n],
	Cdouble d[*n],
	CdoubleCP u,
	doubleCP v,
	CintCP index1,
	CintCP nu)
{
/* Evaluate  v[index1[l]] := spline(u[index1[l]], ...),	    l = 0,..,nu-1
 * Nodes x[index0[i]], coef (y[index0[i]]; b[i],c[i],d[i]); i = 0,..,n-1
 */
	const int n_1 = *n-1;
	int i, j, k, l;
	double ul, dx, tmp;
	i = 0;
	for (l = 0; l < *nu; l++) {
		ul = u[index1[l]];
		if ( ul < x[index0[i]] || (i < n_1 && x[index0[i+1]] < ul) ) {
			/* reset i  such that  x[index0[i]] <= ul <= x[index0[i+1]] : */
			i = 0;
			j = *n;
			do {
				k = (i+j)/2;
				if (ul < x[index0[k]]) j = k;
				else i = k;
			} while (j > i+1);
		}
		dx = ul-x[index0[i]];
		/* for natural splines extrapolate linearly left */
		tmp = (*method == 0 && ul < x[index0[0]]) ? 0.0 : d[i];
		v[index1[l]] = y[index0[i]]+dx*( b[i]+dx*(c[i]+dx*tmp) );
	}
	return;
} // spline_evalI
