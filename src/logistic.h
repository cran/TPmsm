
typedef struct {
	int n, *IPIV;
	double *B, *U;
	int lwork;
	double *F, *W;
} logitW;

logitW	*logitW_Create(CintCP n);
void logitW_Delete(logitW *WORK);
void predict_logit(CintCP len, Cint subset[*len], CintCP Y, doubleCP P, CintCP n, doubleCP X[*n], CintCP maxit, CdoubleCP epsilon, intCP conv, logitW *WORK);

