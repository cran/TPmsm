
#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdlib.h>
#include <R_ext/Visibility.h>
#include <Rinternals.h>
#include <Rversion.h>
#include "defines.h"
#include "BtoTPCmsm.h"
#include "BtoTPmsm.h"
#include "dgpTP.h"
#include "object.h"
#include "RngStream.h"
#include "RngArray.h"
#include "rthreads.h"
#include "rseed.h"
#include "sample.h"
#include "toTPCmsm.h"
#include "toTPmsm.h"
#include "transAJ.h"
#include "transIPCW1.h"
#include "transIPCW2.h"
#include "transKMPW1.h"
#include "transKMPW2.h"
#include "transKMW.h"
#include "transLIN1.h"
#include "transLIN2.h"
#include "transLS.h"
#include "transPAJ.h"
#include "uniqueTP.h"

#define PREFIX Rf_

#define STRINGIFY2(x) #x
#define STRINGIFY(x) STRINGIFY2(x)
#define PASTE2(a, b) a##b
#define PASTE(a, b) PASTE2(a, b)

#define CALLDEF(name, n)  {STRINGIFY(PASTE(PREFIX, name)), (DL_FUNC) &name, n}

static const R_CallMethodDef CallEntries[] = {
	CALLDEF(BtoTPCmsm, 10),
	CALLDEF(BtoTPmsm1222, 8),
	CALLDEF(BtoTPmsm1323, 8),
	CALLDEF(dgpTP, 7),
	CALLDEF(rsample, 1),
	CALLDEF(rset_num_threads, 1),
	CALLDEF(rset_package_seed, 1),
	CALLDEF(rset_seed, 1),
	CALLDEF(SetClass, 2),
	CALLDEF(toTPCmsm, 7),
	CALLDEF(toTPmsm1222, 5),
	CALLDEF(toTPmsm1323, 5),
	CALLDEF(TransPROBAJ, 3),
	CALLDEF(TransPROBIPCW1, 4),
	CALLDEF(TransPROBIPCW2, 8),
	CALLDEF(TransPROBKMPW1, 4),
	CALLDEF(TransPROBKMPW2, 4),
	CALLDEF(TransPROBKMW, 4),
	CALLDEF(TransPROBLIN1, 3),
	CALLDEF(TransPROBLIN2, 7),
	CALLDEF(TransPROBLS, 9),
	CALLDEF(TransPROBPAJ, 3),
	CALLDEF(uniqueTIME, 3),
	CALLDEF(uniqueCOV, 2),
	{NULL, NULL, 0}
};

void attribute_visible R_init_TPmsm(DllInfo *dll) {
	R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
	#if defined(R_VERSION) && R_VERSION >= R_Version(2, 16, 0)
	R_forceSymbols(dll, TRUE);
	#endif

	#ifdef _OPENMP
	global_num_procs = omp_get_num_procs();
	set_num_threads(global_num_procs);
	#endif

	RngArray = (RngStream*)malloc( global_num_procs*sizeof(RngStream) );
	if (RngArray == NULL) error("R_init_TPmsm: No more memory\n");
	RngArray_CreateStream(&global_num_threads, RngArray);

	SEXP TPmsm_NS = R_FindNamespace( mkString("TPmsm") );
	if (TPmsm_NS == R_UnboundValue) error("missing 'TPmsm' namespace: should never happen");
	if ( !isEnvironment(TPmsm_NS) ) error("'TPmsm' namespace not determined correctly");
	return;
} // R_init_TPmsm
