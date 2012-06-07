
#include <stdlib.h>
#include <Rdefines.h>
#include "defines.h"

#define toTPmsm0 \
	SEXP a3d, h; \
	a3d = VECTOR_ELT(lst, 0); \
	h = VECTOR_ELT(lst, 1); \
	register int i; \
	int nt = GET_LENGTH(UT); \
	SEXP mest; \
	PROTECT( mest = allocMatrix(REALSXP, nt, 5) ); \

#define toTPmsm1 \
	const char *name1 = CHAR( STRING_ELT(statenames, 0) ); \
	const char *name2 = CHAR( STRING_ELT(statenames, 1) ); \
	const char *name3 = CHAR( STRING_ELT(statenames, 2) ); \
	int j, k; \
	i = strlen(name1), j = strlen(name2), k = strlen(name3); \
	char *name11 = (char*)malloc( (i+i+2)*sizeof(char) ); \
	strcpy(name11, name1); \
	strcpy(&name11[i], " "); \
	strcpy(&name11[i+1], name1); \
	char *name12 = (char*)malloc( (i+j+2)*sizeof(char) ); \
	strcpy(name12, name1); \
	strcpy(&name12[i], " "); \
	strcpy(&name12[i+1], name2); \
	char *name13 = (char*)malloc( (i+k+2)*sizeof(char) ); \
	strcpy(name13, name1); \
	strcpy(&name13[i], " "); \
	strcpy(&name13[i+1], name3); \
	char *name22 = (char*)malloc( (j+j+2)*sizeof(char) ); \
	strcpy(name22, name2); \
	strcpy(&name22[j], " "); \
	strcpy(&name22[j+1], name2); \
	char *name23 = (char*)malloc( (j+k+2)*sizeof(char) ); \
	strcpy(name23, name2); \
	strcpy(&name23[j], " "); \
	strcpy(&name23[j+1], name3); \
	SEXP colnames; \
	PROTECT( colnames = NEW_CHARACTER(5) ); \
	SET_STRING_ELT( colnames, 0, mkChar(name11) ); \
	SET_STRING_ELT( colnames, 1, mkChar(name12) ); \
	SET_STRING_ELT( colnames, 2, mkChar(name13) ); \
	SET_STRING_ELT( colnames, 3, mkChar(name22) ); \
	SET_STRING_ELT( colnames, 4, mkChar(name23) ); \
	free(name11); free(name12); free(name13); \
	free(name22); free(name23); \
	SEXP dimnames; \
	PROTECT( dimnames = NEW_LIST(2) ); \
	SET_ELEMENT(dimnames, 0, R_NilValue); \
	SET_ELEMENT(dimnames, 1, colnames); \
	SET_DIMNAMES(mest, dimnames); \
	SEXP list; \
	PROTECT( list = NEW_LIST(11) ); \
	SET_ELEMENT( list, 0, GET_CLASS(lst) ); \
	SET_ELEMENT(list, 1, mest); \
	SET_ELEMENT(list, 2, R_NilValue); \
	SET_ELEMENT(list, 3, R_NilValue); \
	SET_ELEMENT(list, 4, UT); \
	SET_ELEMENT(list, 5, s); \
	SET_ELEMENT(list, 6, t); \
	SET_ELEMENT(list, 7, h); \
	SET_ELEMENT(list, 8, statenames); \
	SET_ELEMENT(list, 9, R_NilValue); \
	SET_ELEMENT(list, 10, R_NilValue); \
	SEXP names; \
	PROTECT( names = NEW_CHARACTER(11) ); \
	SET_STRING_ELT( names, 0, mkChar("method") ); \
	SET_STRING_ELT( names, 1, mkChar("est") ); \
	SET_STRING_ELT( names, 2, mkChar("inf") ); \
	SET_STRING_ELT( names, 3, mkChar("sup") ); \
	SET_STRING_ELT( names, 4, mkChar("time") ); \
	SET_STRING_ELT( names, 5, mkChar("s") ); \
	SET_STRING_ELT( names, 6, mkChar("t") ); \
	SET_STRING_ELT( names, 7, mkChar("h") ); \
	SET_STRING_ELT( names, 8, mkChar("state.names") ); \
	SET_STRING_ELT( names, 9, mkChar("n.boot") ); \
	SET_STRING_ELT( names, 10, mkChar("conf.level") ); \
	SET_NAMES(list, names); \
	SEXP classl; \
	PROTECT( classl = NEW_CHARACTER(1) ); \
	SET_STRING_ELT( classl, 0, mkChar("TPmsm") ); \
	SET_CLASS(list, classl); \
	UNPROTECT(6); \
	return list; \

SEXP toTPmsm1222(
	SEXP lst,
	SEXP UT,
	SEXP s,
	SEXP t,
	SEXP statenames)
{
	toTPmsm0
	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for (i = 0; i < nt; i++) {
		if (REAL(a3d)[i] < 0) REAL(mest)[i] = 0;
		else if (REAL(a3d)[i] > 1) REAL(mest)[i] = 1;
		else REAL(mest)[i] = REAL(a3d)[i];
		if (REAL(a3d)[i+nt] < 0) REAL(mest)[i+nt] = 0;
		else if (REAL(a3d)[i+nt] > 1) REAL(mest)[i+nt] = 1;
		else REAL(mest)[i+nt] = REAL(a3d)[i+nt];
		REAL(mest)[i+nt*2] = 1-REAL(mest)[i]-REAL(mest)[i+nt];
		if (REAL(mest)[i+nt*2] < 0) {
			REAL(mest)[i+nt] = 1-REAL(mest)[i];
			REAL(mest)[i+nt*2] = 0;
		}
		if (REAL(a3d)[i+nt*2] < 0) REAL(mest)[i+nt*3] = 0;
		else if (REAL(a3d)[i+nt*2] > 1) REAL(mest)[i+nt*3] = 1;
		else REAL(mest)[i+nt*3] = REAL(a3d)[i+nt*2];
		REAL(mest)[i+nt*4] = 1-REAL(mest)[i+nt*3];
	}
	toTPmsm1
} // toTPmsm1222

SEXP toTPmsm1323(
	SEXP lst,
	SEXP UT,
	SEXP s,
	SEXP t,
	SEXP statenames)
{
	toTPmsm0
	#ifdef _OPENMP
	#pragma omp parallel for private(i)
	#endif
	for (i = 0; i < nt; i++) {
		if (REAL(a3d)[i] < 0) REAL(mest)[i] = 0;
		else if (REAL(a3d)[i] > 1) REAL(mest)[i] = 1;
		else REAL(mest)[i] = REAL(a3d)[i];
		if (REAL(a3d)[i+nt] < 0) REAL(mest)[i+nt*2] = 0;
		else if (REAL(a3d)[i+nt] > 1) REAL(mest)[i+nt*2] = 1;
		else REAL(mest)[i+nt*2] = REAL(a3d)[i+nt];
		REAL(mest)[i+nt] = 1-REAL(mest)[i]-REAL(mest)[i+nt*2];
		if (REAL(mest)[i+nt] < 0) {
			REAL(mest)[i+nt*2] = 1-REAL(mest)[i];
			REAL(mest)[i+nt] = 0;
		}
		if (REAL(a3d)[i+nt*2] < 0) REAL(mest)[i+nt*4] = 0;
		else if (REAL(a3d)[i+nt*2] > 1) REAL(mest)[i+nt*4] = 1;
		else REAL(mest)[i+nt*4] = REAL(a3d)[i+nt*2];
		REAL(mest)[i+nt*3] = 1-REAL(mest)[i+nt*4];
	}
	toTPmsm1
} // toTPmsm1323
