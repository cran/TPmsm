
#include <Rdefines.h>
#include "defines.h" // needed by "RngArray.h"
#include "RngStream.h"
#include "RngArray.h"

int global_num_procs = 1;
int global_num_threads = 1;

void set_num_threads(int tnum) {
	global_num_threads = tnum > global_num_procs ? global_num_procs : tnum;
	return;
} // set_num_threads

SEXP rset_num_threads(SEXP arg_snum) {
	SEXP ret_snum;
	PROTECT( ret_snum = NEW_INTEGER(1) );
	*INTEGER(ret_snum) = global_num_threads;
	if ( !isNull(arg_snum) ) {
		set_num_threads( *INTEGER(arg_snum) );
		RngArray_DeleteStream(&global_num_procs, RngArray);
		RngArray_CreateStream(&global_num_threads, RngArray);
	}
	UNPROTECT(1);
	return ret_snum;
} // rset_num_threads
