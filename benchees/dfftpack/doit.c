/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "dfftpack")
BENCH_DOC("author", "Paul N. Swarztrauber and Hugh C. Pumphrey")
BENCH_DOC("year", "1985")
BENCH_DOC("version", "1.0")
BENCH_DOC("language", "FORTRAN")
END_BENCH_DOC


int can_do(struct problem *p)
{
     return (sizeof(bench_real) == sizeof(double) && p->rank == 1);
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_1d_fftpack(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_fftpack(p, in, -1.0);
}

static bench_real *WSAVE;

extern void F77_FUNC(zffti, ZFFTI)();
extern void F77_FUNC(dffti, DFFTI)();
extern void F77_FUNC(zfftf, ZFFTF)();
extern void F77_FUNC(dfftf, DFFTF)();
extern void F77_FUNC(zfftb, ZFFTB)();
extern void F77_FUNC(dfftb, DFFTB)();

void setup(struct problem *p)
{
     int n;
 
     BENCH_ASSERT(can_do(p));
     n = p->n[0];
 
     if (p->kind == PROBLEM_COMPLEX) {
	  WSAVE = bench_malloc((4 * n + 15) * sizeof(bench_real));
	  F77_FUNC(zffti, ZFFTI)(&n, WSAVE);
     } else {
	  WSAVE = bench_malloc((2 * n + 15) * sizeof(bench_real));
	  F77_FUNC(dffti, DFFTI)(&n, WSAVE);
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     int n = p->n[0];
     void *in = p->in;
     bench_real *wsave = WSAVE;

     if (p->kind == PROBLEM_COMPLEX) {
	  if (p->sign == -1) {
	       for (i = 0; i < iter; ++i) {
		    F77_FUNC(zfftf, ZFFTF)(&n, in, wsave);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) {
		    F77_FUNC(zfftb, ZFFTB)(&n, in, wsave);
	       }
	  }
     } else {
	  if (p->sign == -1) {
	       for (i = 0; i < iter; ++i) {
		    F77_FUNC(dfftf, DFFTF)(&n, in, wsave);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) {
		    F77_FUNC(dfftb, DFFTB)(&n, in, wsave);
	       }
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(WSAVE);
}
