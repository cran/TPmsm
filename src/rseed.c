
#include <Rdefines.h>
#include "defines.h" // needed by "RngArray.h"
#include "RngStream.h"
#include "RngArray.h"
#include "rthreads.h"

SEXP rset_package_seed(SEXP arg_seed)
{
	register int i;
	unsigned long seed[6];
	for (i = 0; i < 6; i++) {
		seed[i] = (unsigned long)REAL(arg_seed)[i];
	}
	RngStream_SetPackageSeed(seed);
	RngArray_DeleteStream(&global_num_procs, RngArray);
	RngArray_CreateStream(&global_num_threads, RngArray);
	return R_NilValue;
} // rset_package_seed

SEXP rset_seed(SEXP arg_seed)
{
	register int i, j;
	unsigned long seed[6];
	SEXP ret_seed, rng_seed[global_num_threads], classl;
	PROTECT( ret_seed = NEW_LIST(global_num_threads) );
	for (i = 0; i < global_num_threads; i++) {
		PROTECT( rng_seed[i] = NEW_NUMERIC(6) );
		RngStream_GetState(RngArray[i], seed);
		for (j = 0; j < 6; j++) REAL(rng_seed[i])[j] = seed[j];
		SET_ELEMENT(ret_seed, i, rng_seed[i]);
	}
	if ( !isNull(arg_seed) ) {
		for (i = 0; i < global_num_threads; i++) {
			for (j = 0; j < 6; j++)
				seed[j]= (unsigned long)REAL( VECTOR_ELT(arg_seed, i) )[j];
			RngStream_SetSeed(RngArray[i], seed);
		}
	}
	PROTECT( classl = NEW_CHARACTER(1) );
	SET_STRING_ELT( classl, 0, mkChar("TPmsmSeed") );
	SET_CLASS(ret_seed, classl);
	UNPROTECT(global_num_threads+2);
	return ret_seed;
} // rset_seed
