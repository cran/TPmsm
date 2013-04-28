
typedef const double *const CdoubleCP;
typedef double *const doubleCP;
typedef const double Cdouble;
typedef const int *const CintCP;
typedef int *const intCP;
typedef const int Cint;
typedef const char *const CcharCP;

typedef enum {
	SINT_PTR, // short int pointer
	INT_PTR, // int pointer
	REAL_PTR, // double pointer
} etype;

typedef union {
	short int *shortinteger;
	int *integer;
	double *real;
} utype;

typedef struct {
	etype type; // type of pointer
	int length; // length of array
	utype ptr; // pointer to array
} stype;

typedef const stype *const CstypeCP;
typedef stype *const stypeCP;
