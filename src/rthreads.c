
#include <Rdefines.h>
#ifdef _OPENMP
#include <omp.h>
#endif

int global_num_threads = 1;

void set_num_threads(int tnum) {
	#ifdef _OPENMP
	int pnum = omp_get_num_procs();
	#else
	int pnum = 1;
	#endif
	global_num_threads = tnum > pnum ? pnum : tnum;
	return;
} // set_num_threads

SEXP rset_num_threads(SEXP arg_snum) {
	SEXP ret_snum;
	PROTECT( ret_snum = NEW_INTEGER(1) );
	*INTEGER(ret_snum) = global_num_threads;
	if ( !isNull(arg_snum) ) set_num_threads( *INTEGER(arg_snum) );
	UNPROTECT(1);
	return ret_snum;
} // rset_num_threads
