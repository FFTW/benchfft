#include "bench-user.h"

#define NAME "nag-multiple"
#define NOTE "Uses functions c06frf and friends: interface for multiple simultaneous transforms (not used), precomputed trig. tables, and extra workspace supplied."

#define ONLY_1D 1

#define C06FRF F77_FUNC(c06frf,C06FRF)
extern void C06FRF(unsigned int *m, unsigned int *n, double *x, double *y, char *init,
		   double *trig, double *work, int *ifail);
#define fft(n,x,y,trig,work) C06FRF(&one, &n, x, y, &SUBS, trig, work, &ifail)
#define init(n,x,y,trig,work) C06FRF(&one, &n, x, y, &INIT, trig, work, &ifail)

#define C06FPF F77_FUNC(c06fpf,C06FPF)
extern void C06FPF(unsigned int *m, unsigned int *n, double *x, char *init,
		   double *trig, double *work, int *ifail);
#define fftrc(n,x,t,work) C06FPF(&one, &n, x, &SUBS, t, work, &ifail)
#define initrc(n,x,t,work) C06FPF(&one, &n, x, &INIT, t, work, &ifail)

#define C06FQF F77_FUNC(c06fqf,C06FQF)
extern void C06FQF(unsigned int *m, unsigned int *n, double *x, char *init,
		   double *trig, double *work, int *ifail);
#define fftcr(n,x,t,work) C06FQF(&one, &n, x, &SUBS, t, work, &ifail)
#define initcr(n,x,t,work) C06FQF(&one, &n, x, &INIT, t, work, &ifail)

#include "doit.c"
