#include "bench-user.h"

#define NAME "nag-inplace"
#define NOTE "Uses c06ecf and friends (in-place, no workspace)."

#define ONLY_1D 1

#define C06ECF F77_FUNC(c06ecf,C06ECF)
extern void C06ECF(double *x, double *y, unsigned int *n, int *ifail);
#define fft(n,x,y,trig,work) C06ECF(x, y, &n, &ifail)
#define init(n,x,y,trig,work) 0

#define C06EAF F77_FUNC(c06eaf,C06EAF)
extern void C06EAF(double *x, unsigned int *n, int *ifail);
#define fftrc(n,x,trig,work) C06EAF(x, &n, &ifail)
#define initrc(n,x,trig,work) 0

#define C06EBF F77_FUNC(c06ebf,C06EBF)
extern void C06EBF(double *x, unsigned int *n, int *ifail);
#define fftcr(n,x,trig,work) C06EBF(x, &n, &ifail)
#define initcr(n,x,trig,work) 0

#include "doit.c"
