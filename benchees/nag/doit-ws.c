#include "bench-user.h"

#define NAME "nag-workspace"
#define NOTE "Uses c06fcf and friends (extra workspace supplied)."

#define C06FCF F77_FUNC(c06fcf,C06FCF)
extern void C06FCF(double *x, double *y, unsigned int *n,
		   double *work, int *ifail);
#define fft(n,x,y,trig,work) C06FCF(x, y, &n, work, &ifail)
#define init(n,x,y,trig,work) 0

#define C06FAF F77_FUNC(c06faf,C06FAF)
extern void C06FAF(double *x, unsigned int *n, double *work, int *ifail);
#define fftrc(n,x,trig,work) C06FAF(x, &n, work, &ifail)
#define initrc(n,x,trig,work) 0

#define C06FBF F77_FUNC(c06fbf,C06FBF)
extern void C06FBF(double *x, unsigned int *n, double *work, int *ifail);
#define fftcr(n,x,trig,work) C06FBF(x, &n, work, &ifail)
#define initcr(n,x,trig,work) 0

#include "doit.c"
