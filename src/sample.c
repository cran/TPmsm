
#ifdef _OPENMP
#include <omp.h>
#endif
#include <Rdefines.h>
#include "defines.h" // needed by "RngArray.h"
#include "RngStream.h"
#include "RngArray.h"
#include "rthreads.h"

SEXP rsample(SEXP arg_n)
{
	SEXP ret_i;
	PROTECT( ret_i = NEW_INTEGER( *INTEGER(arg_n) ) );
	#ifdef _OPENMP
	#pragma omp parallel num_threads(global_num_threads)
	#endif
	{
		register int i;
		int t;
		#ifdef _OPENMP
		t = omp_get_thread_num();
		#else
		t = 0;
		#endif
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for (i = 0; i < *INTEGER(arg_n); i++) {
			INTEGER(ret_i)[i] = (int)RngStream_RandInt( RngArray[t], 1, (long)*INTEGER(arg_n) );
		}
	}
	UNPROTECT(1);
	return ret_i;
} // rsample
