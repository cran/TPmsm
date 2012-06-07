
void (*kchar2ptr(SEXP window))(CdoubleCP, CintCP, CdoubleCP, CdoubleCP, doubleCP);
void NWWeights( CdoubleCP X, CintCP n, CdoubleCP x, CdoubleCP h, doubleCP W, void (*kfunc)(CdoubleCP, CintCP, CdoubleCP, CdoubleCP, doubleCP) );
void LLWeights( CdoubleCP X, CintCP n, CdoubleCP x, CdoubleCP h, doubleCP W, void (*kfunc)(CdoubleCP, CintCP, CdoubleCP, CdoubleCP, doubleCP) );
