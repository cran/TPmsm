
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <Rdefines.h>
#include "defines.h"
#include "quantile.h"
#include "rthreads.h"

SEXP BtoTPCmsm(
	SEXP lst,
	SEXP UT,
	SEXP UX,
	SEXP s,
	SEXP t,
	SEXP x,
	SEXP statenames,
	SEXP nboot,
	SEXP conflevel,
	SEXP methodboot)
{
	SEXP a4d, h;
	a4d = VECTOR_ELT(lst, 0);
	h = VECTOR_ELT(lst, 1);
	register int i, j, k;
	register int64_t y;
	const int nt = GET_LENGTH(UT), nx = GET_LENGTH(UX);
	const int64_t ntx = nt*nx;
	SEXP aest;
	PROTECT( aest = alloc3DArray(REALSXP, nt, nx, 5) );
	int n = 2, nb = INTEGER( GET_DIM(a4d) )[0];
	double P[n], Q[n];
	P[0] = ( 1-*REAL(conflevel) )/2;
	P[1] = ( 1+*REAL(conflevel) )/2;
	SEXP ainf, asup;
	PROTECT( ainf = alloc3DArray(REALSXP, nt, nx, 5) );
	PROTECT( asup = alloc3DArray(REALSXP, nt, nx, 5) );
	#ifdef _OPENMP
	#pragma omp parallel for num_threads(global_num_threads) private(i, j, k, y, Q)
	#endif
	for (i = 0; i < nt; i++) {
		for (y = i, j = 0; j < nx; j++) {
			for (k = 0; k < 4; k++) {
				REAL(aest)[y] = REAL(a4d)[nb*y]; // must be computed before quantile_d
				quantile_d(&nb, &REAL(a4d)[nb*y], &n, P, Q); // sorts the array internally
				REAL(ainf)[y] = Q[0];
				REAL(asup)[y] = Q[1];
				y += ntx; // y = i+nt*j+ntx*k
			}
			REAL(aest)[y] = 1-REAL(aest)[y-ntx]; // y = i+nt*j+ntx*4, y-ntx = i+nt*j+ntx*3
			REAL(ainf)[y] = 1-Q[1];
			REAL(asup)[y] = 1-Q[0];
			y -= ntx*4; // y = i+nt*j+ntx*4-ntx*4 = i+nt*j+ntx*0 = i+nt*j
			y += nt; // y = i+nt*j
		}
	}
	if (strcmp(CHAR( STRING_ELT(methodboot, 0) ), "basic") == 0) {
		#ifdef _OPENMP
		#pragma omp parallel for num_threads(global_num_threads) private(i, j, k, y, P)
		#endif
		for (i = 0; i < nt; i++) {
			for (y = i, j = 0; j < nx; j++) {
				for (k = 0; k < 5; k++) {
					P[0] = 2*REAL(aest)[y]-REAL(asup)[y];
					if (P[0] < 0) REAL(ainf)[y] = 0;
					else if (P[0] > 1) REAL(ainf)[y] = 1;
					else REAL(ainf)[y] = P[0];
					P[1] = 2*REAL(aest)[y]-REAL(ainf)[y];
					if (P[1] < 0) REAL(asup)[y] = 0;
					else if (P[1] > 1) REAL(asup)[y] = 1;
					else REAL(asup)[y] = P[1];
					y += ntx; // y = i+nt*j+ntx*k
				}
				y -= ntx*5; // y = i+nt*j+ntx*5-ntx*5 = i+nt*j+ntx*0 = i+nt*j
				y += nt; // y = i+nt*j
			}
		}
	}
	const char *name1 = CHAR( STRING_ELT(statenames, 0) );
	const char *name2 = CHAR( STRING_ELT(statenames, 1) );
	const char *name3 = CHAR( STRING_ELT(statenames, 2) );
	i = strlen(name1), j = strlen(name2), k = strlen(name3);
	char *name11 = (char*)malloc( (i+i+2)*sizeof(char) );
	if (name11 == NULL) error("BtoTPCmsm: No more memory\n");
	strcpy(name11, name1);
	strcpy(&name11[i], " ");
	strcpy(&name11[i+1], name1);
	char *name12 = (char*)malloc( (i+j+2)*sizeof(char) );
	if (name12 == NULL) error("BtoTPCmsm: No more memory\n");
	strcpy(name12, name1);
	strcpy(&name12[i], " ");
	strcpy(&name12[i+1], name2);
	char *name13 = (char*)malloc( (i+k+2)*sizeof(char) );
	if (name13 == NULL) error("BtoTPCmsm: No more memory\n");
	strcpy(name13, name1);
	strcpy(&name13[i], " ");
	strcpy(&name13[i+1], name3);
	char *name22 = (char*)malloc( (j+j+2)*sizeof(char) );
	if (name22 == NULL) error("BtoTPCmsm: No more memory\n");
	strcpy(name22, name2);
	strcpy(&name22[j], " ");
	strcpy(&name22[j+1], name2);
	char *name23 = (char*)malloc( (j+k+2)*sizeof(char) );
	if (name23 == NULL) error("BtoTPCmsm: No more memory\n");
	strcpy(name23, name2);
	strcpy(&name23[j], " ");
	strcpy(&name23[j+1], name3);
	SEXP facenames;
	PROTECT( facenames = NEW_CHARACTER(5) );
	SET_STRING_ELT( facenames, 0, mkChar(name11) );
	SET_STRING_ELT( facenames, 1, mkChar(name12) );
	SET_STRING_ELT( facenames, 2, mkChar(name13) );
	SET_STRING_ELT( facenames, 3, mkChar(name22) );
	SET_STRING_ELT( facenames, 4, mkChar(name23) );
	free(name11); free(name12); free(name13);
	free(name22); free(name23);
	SEXP dimnames;
	PROTECT( dimnames = NEW_LIST(3) );
	SET_ELEMENT(dimnames, 0, R_NilValue);
	SET_ELEMENT(dimnames, 1, R_NilValue);
	SET_ELEMENT(dimnames, 2, facenames);
	SET_DIMNAMES(aest, dimnames);
	SET_DIMNAMES(ainf, dimnames);
	SET_DIMNAMES(asup, dimnames);
	SEXP list;
	PROTECT( list = NEW_LIST(13) );
	SET_ELEMENT( list, 0, GET_CLASS(lst) );
	SET_ELEMENT(list, 1, aest);
	SET_ELEMENT(list, 2, ainf);
	SET_ELEMENT(list, 3, asup);
	SET_ELEMENT(list, 4, UT);
	SET_ELEMENT(list, 5, UX);
	SET_ELEMENT(list, 6, s);
	SET_ELEMENT(list, 7, t);
	SET_ELEMENT(list, 8, x);
	SET_ELEMENT(list, 9, h);
	SET_ELEMENT(list, 10, statenames);
	SET_ELEMENT(list, 11, nboot);
	SET_ELEMENT(list, 12, conflevel);
	SEXP names;
	PROTECT( names = NEW_CHARACTER(13) );
	SET_STRING_ELT( names, 0, mkChar("method") );
	SET_STRING_ELT( names, 1, mkChar("est") );
	SET_STRING_ELT( names, 2, mkChar("inf") );
	SET_STRING_ELT( names, 3, mkChar("sup") );
	SET_STRING_ELT( names, 4, mkChar("time") );
	SET_STRING_ELT( names, 5, mkChar("covariate") );
	SET_STRING_ELT( names, 6, mkChar("s") );
	SET_STRING_ELT( names, 7, mkChar("t") );
	SET_STRING_ELT( names, 8, mkChar("x") );
	SET_STRING_ELT( names, 9, mkChar("h") );
	SET_STRING_ELT( names, 10, mkChar("state.names") );
	SET_STRING_ELT( names, 11, mkChar("n.boot") );
	SET_STRING_ELT( names, 12, mkChar("conf.level") );
	SET_NAMES(list, names);
	SEXP classl;
	PROTECT( classl = NEW_CHARACTER(1) );
	SET_STRING_ELT( classl, 0, mkChar("TPCmsm") );
	SET_CLASS(list, classl);
	UNPROTECT(8);
	return list;
} // BtoTPCmsm
