
#include <R_ext/Random.h>
#include <Rdefines.h>
#include <Rmath.h>
#include "defines.h"

typedef void (*Tfunc)(CdoubleCP, CdoubleCP, doubleCP, doubleCP);
typedef void (*Cfunc)(CdoubleCP, doubleCP);

static void expt(
	CdoubleCP pcorr,
	CdoubleCP pdistpar,
	doubleCP t1,
	doubleCP t2)
{
	double u1, u2, v, a;
	u1 = runif(0, 1);
	v = runif(0, 1);
	a = *pcorr*(2*u1-1);
	u2 = 2*v/( 1-a+sqrt(R_pow_di(1-a, 2)+4*a*v) );
	*t1 = -pdistpar[0]*log(1-u1);
	*t2 = -pdistpar[1]*log(1-u2);
	return;
} // expt

static void weibullt(
	CdoubleCP pcorr,
	CdoubleCP pdistpar,
	doubleCP t1,
	doubleCP t2)
{
	register int i;
	double u[5], v;
	for (i = 0; i < 5; i++) u[i] = runif(0, 1);
	if (u[4] > *pcorr) v = -log(u[3]);
	else v = -log(u[1])-log(u[2]);
	*t1 = R_pow(u[0], *pcorr/pdistpar[0])*R_pow(v, 1/pdistpar[0])*pdistpar[1];
	*t2 = R_pow(1-u[0], *pcorr/pdistpar[2])*R_pow(v, 1/pdistpar[2])*pdistpar[3];
	return;
} // weibullt

static void runif0(
	CdoubleCP pcenspar,
	doubleCP c)
{
	*c = runif(0, *pcenspar);
	return;
} // runif0

static void rexp0(
	CdoubleCP pcenspar,
	doubleCP c)
{
	double u;
	u = runif(0, 1);
	*c = -*pcenspar*log(1-u);
	return;
} // rexp0

static void cens2(
	doubleCP pT1,
	intCP pE1,
	doubleCP pS,
	intCP pE,
	CintCP pn,
	Tfunc tfunc,
	CdoubleCP pcorr,
	CdoubleCP pdistpar,
	Cfunc cfunc,
	CdoubleCP pcenspar,
	CdoubleCP pstate2prob)
{
	register int i;
	double t1, t2, c, b;
	for (i = 0; i < *pn; i++) {
		cfunc(pcenspar, &c);
		tfunc(pcorr, pdistpar, &t1, &t2);
		b = rbinom(1, *pstate2prob);
		pT1[i] = fmin2(t1, c);
		pE1[i] = (t1 <= c);
		pS[i] = pT1[i]+b*pE1[i]*fmin2(t2, c-t1);
		pE[i] = (1-b)*pE1[i]+b*(t2 <= c-t1);
	}
	return;
} // cens2

/*
Author:
	Artur Agostinho Araujo <artur.stat@gmail.com>

Description:
	Generates bivariate censored gap times from some known copula functions.

Parameters:
	n				sample size.
	corr			correlation parameter.
	dist			distribution.
	distpar			vector of parameters for the distribution.
	modelcens		model for censorship.
	censpar			parameter for the censorship distribution.
	state2prob		the proportion of individuals entering state2.

Return value:
	Returns an object of class 'survTP'.
*/

SEXP dgpTP(
	SEXP n,
	SEXP corr,
	SEXP dist,
	SEXP distpar,
	SEXP modelcens,
	SEXP censpar,
	SEXP state2prob)
{
	CintCP pn = INTEGER_POINTER(n);
	CdoubleCP pcorr = NUMERIC_POINTER(corr);
	CcharCP pdist = CHAR( STRING_ELT(dist, 0) );
	CdoubleCP pdistpar = NUMERIC_POINTER(distpar);
	CcharCP pmodelcens = CHAR( STRING_ELT(modelcens, 0) );
	CdoubleCP pcenspar = NUMERIC_POINTER(censpar);
	CdoubleCP pstate2prob = NUMERIC_POINTER(state2prob);
	if (*pn <= 0) error("Argument 'n' must be greater than zero");
	CcharCP model1 = "uniform", model2 = "exponential";
	CcharCP model3 = "weibull";
	if ( !(strcmp(pdist, model2) == 0 || strcmp(pdist, model3) == 0) ) error("Argument 'dist' must be one of 'weibull' or 'exponential'");
	Tfunc tfunc = expt;
	if (strcmp(pdist, model2) == 0) {
		if (*pcorr < -1 || *pcorr > 1) error("Argument 'corr' with dist='exponential' must be greater or equal to -1 and lower or equal to 1");
		if (GET_LENGTH(distpar) != 2) error("Argument 'dist.par' with 'dist=exponential' must be a vector with lenght 2");
		if (pdistpar[0] <= 0 || pdistpar[1] <= 0) error("Argument 'dist.par' must be greater than 0");
	}
	if (strcmp(pdist, model3) == 0) {
		if (*pcorr <= 0 || *pcorr > 1) error("Argument 'corr' with 'dist=weibull' must be greater than 0 and lower or equal to 1");
		if (GET_LENGTH(distpar) != 4) error("Argument 'dist.par' with 'dist=weibull' must be a vector with lenght 4");
		if (pdistpar[0] <= 0 || pdistpar[1] <= 0 || pdistpar[2] <= 0 || pdistpar[3] <= 0) error("Argument 'dist.par' must be greater than 0");
		tfunc = weibullt;
	}
	if ( !(strcmp(pmodelcens, model1) == 0 || strcmp(pmodelcens, model2) == 0) ) error("Argument 'model.cens' must be one of 'uniform' or 'exponential'");
	Cfunc cfunc = runif0;
	if (strcmp(pmodelcens, model1) == 0) {
		if (*pcenspar < 0) error("Argument 'cens.par' with 'model.cens=uniform' must be greater or equal than 0");
	}
	if (strcmp(pmodelcens, model2) == 0) {
		if (*pcenspar <= 0) error("Argument 'cens.par' with 'model.cens=exponential' must be greater than 0");
		cfunc = rexp0;
	}
	if (*pstate2prob < 0 || *pstate2prob > 1) error("Argument 'state2.prob' must be greater or equal to 0 and lower or equal to 1");
	void (*func)(doubleCP, intCP, doubleCP, intCP, CintCP, Tfunc, CdoubleCP, CdoubleCP, Cfunc, CdoubleCP, CdoubleCP);
	func = cens2;
	SEXP T1, E1, S, E;
	PROTECT( T1 = NEW_NUMERIC(*pn) );
	PROTECT( E1 = NEW_INTEGER(*pn) );
	PROTECT( S = NEW_NUMERIC(*pn) );
	PROTECT( E = NEW_INTEGER(*pn) );
	GetRNGstate();
	func(REAL(T1), INTEGER(E1), REAL(S), INTEGER(E), pn, tfunc, pcorr, pdistpar, cfunc, pcenspar, pstate2prob);
	PutRNGstate();
	SEXP data;
	PROTECT( data = NEW_LIST(4) );
	SET_ELEMENT(data, 0, T1);
	SET_ELEMENT(data, 1, E1);
	SET_ELEMENT(data, 2, S);
	SET_ELEMENT(data, 3, E);
	SEXP names;
	PROTECT( names = NEW_CHARACTER(4) );
	SET_STRING_ELT( names, 0, mkChar("time1") );
	SET_STRING_ELT( names, 1, mkChar("event1") );
	SET_STRING_ELT( names, 2, mkChar("Stime") );
	SET_STRING_ELT( names, 3, mkChar("event") );
	SET_NAMES(data, names);
	SEXP rows;
	PROTECT( rows = NEW_INTEGER(*pn) );
	register int i;
	for (i = 0; i < *pn; i++) INTEGER(rows)[i] = i+1;
	SET_ATTR(data, R_RowNamesSymbol, rows);
	SEXP classdf;
	PROTECT( classdf = NEW_CHARACTER(1) );
	SET_STRING_ELT( classdf, 0, mkChar("data.frame") );
	SET_CLASS(data, classdf);
	SEXP list;
	PROTECT( list = NEW_LIST(1) );
	SET_ELEMENT(list, 0, data);
	SEXP classl;
	PROTECT( classl = NEW_CHARACTER(1) );
	SET_STRING_ELT( classl, 0, mkChar("survTP") );
	SET_CLASS(list, classl);
	UNPROTECT(10);
	return list;
} // dgpTP
