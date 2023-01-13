
#include <R_ext/Visibility.h>
#include <Rdefines.h>
#include <stdlib.h>
#include "defines.h"
#include "RngStream.h"
#include "RngArray.h"
#include "rthreads.h"

void attribute_visible R_unload_TPmsm(void) {
	RngArray_DeleteStream(&global_num_procs, RngArray);
	free(RngArray);
	return;
} // R_unload_TPmsm
