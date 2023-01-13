
#include <stdint.h>
#include <stdlib.h>
#include <Rdefines.h>
#include "defines.h"
#include "rthreads.h"

#define toTPmsm0 \
	SEXP a3d, h; \
	a3d = VECTOR_ELT(lst, 0); \
	h = VECTOR_ELT(lst, 1); \
	register int i, j; \
	register int64_t k; \
	Cint nt = GET_LENGTH(UT); \
	SEXP mest; \
	PROTECT( mest = allocMatrix(REALSXP, nt, 5) ); \

#define toTPmsm1 \
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
	#pragma omp parallel for num_threads(global_num_threads) private(i, j, k)
	#endif
	for (i = 0; i < nt; i++) {
		for (k = i, j = 0; j < 4; j++) {
			REAL(mest)[k] = REAL(a3d)[k];
			k += nt; // k = i+nt*j
		}
		REAL(mest)[k] = 1-REAL(mest)[k-nt]; // k = i+nt*4, k-nt = i+nt*3
	}
	const char *name1 = CHAR( STRING_ELT(statenames, 0) );
	const char *name2 = CHAR( STRING_ELT(statenames, 1) );
	const char *name3 = CHAR( STRING_ELT(statenames, 2) );
	size_t u = strlen(name1), v = strlen(name2), w = strlen(name3);
	char *name11 = (char*)malloc( (u+u+2)*sizeof(char) );
	if (name11 == NULL) error("toTPmsm1222: No more memory\n");
	strcpy(name11, name1);
	strcpy(&name11[u], " ");
	strcpy(&name11[u+1], name1);
	char *name12 = (char*)malloc( (u+v+2)*sizeof(char) );
	if (name12 == NULL) error("toTPmsm1222: No more memory\n");
	strcpy(name12, name1);
	strcpy(&name12[u], " ");
	strcpy(&name12[u+1], name2);
	char *name13 = (char*)malloc( (u+w+2)*sizeof(char) );
	if (name13 == NULL) error("toTPmsm1222: No more memory\n");
	strcpy(name13, name1);
	strcpy(&name13[u], " ");
	strcpy(&name13[u+1], name3);
	char *name22 = (char*)malloc( (v+v+2)*sizeof(char) );
	if (name22 == NULL) error("toTPmsm1222: No more memory\n");
	strcpy(name22, name2);
	strcpy(&name22[v], " ");
	strcpy(&name22[v+1], name2);
	char *name23 = (char*)malloc( (v+w+2)*sizeof(char) );
	if (name23 == NULL) error("toTPmsm1222: No more memory\n");
	strcpy(name23, name2);
	strcpy(&name23[v], " ");
	strcpy(&name23[v+1], name3);
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
	#pragma omp parallel for num_threads(global_num_threads) private(i, j, k)
	#endif
	for (i = 0; i < nt; i++) {
		for (k = i, j = 0; j < 3; j++) {
			REAL(mest)[k] = REAL(a3d)[k];
			k += nt; // k = i+nt*j
		}
		REAL(mest)[k+nt] = REAL(a3d)[k]; // k+nt = i+nt*4, k = i+nt*3
		REAL(mest)[k] = 1-REAL(a3d)[k];
	}
	const char *name1 = CHAR( STRING_ELT(statenames, 0) );
	const char *name2 = CHAR( STRING_ELT(statenames, 1) );
	const char *name3 = CHAR( STRING_ELT(statenames, 2) );
	size_t u = strlen(name1), v = strlen(name2), w = strlen(name3);
	char *name11 = (char*)malloc( (u+u+2)*sizeof(char) );
	if (name11 == NULL) error("toTPmsm1323: No more memory\n");
	strcpy(name11, name1);
	strcpy(&name11[u], " ");
	strcpy(&name11[u+1], name1);
	char *name12 = (char*)malloc( (u+v+2)*sizeof(char) );
	if (name12 == NULL) error("toTPmsm1323: No more memory\n");
	strcpy(name12, name1);
	strcpy(&name12[u], " ");
	strcpy(&name12[u+1], name2);
	char *name13 = (char*)malloc( (u+w+2)*sizeof(char) );
	if (name13 == NULL) error("toTPmsm1323: No more memory\n");
	strcpy(name13, name1);
	strcpy(&name13[u], " ");
	strcpy(&name13[u+1], name3);
	char *name22 = (char*)malloc( (v+v+2)*sizeof(char) );
	if (name22 == NULL) error("toTPmsm1323: No more memory\n");
	strcpy(name22, name2);
	strcpy(&name22[v], " ");
	strcpy(&name22[v+1], name2);
	char *name23 = (char*)malloc( (v+w+2)*sizeof(char) );
	if (name23 == NULL) error("toTPmsm1323: No more memory\n");
	strcpy(name23, name2);
	strcpy(&name23[v], " ");
	strcpy(&name23[v+1], name3);
	toTPmsm1
} // toTPmsm1323
