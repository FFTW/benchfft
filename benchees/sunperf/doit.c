/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>
#include <stdio.h>

#ifdef HAVE_SUNPERF_H
#include <sunperf.h>
#endif

static const char *mkvers(void)
{
#ifdef HAVE_SUNPERF_VERSION_
     int a, b, c;
     static char *b[128];

     /* 
	using FORTRAN interface.  The documented C interface cannot
	possibly work because it takes a, b, c by value and it returns
	void
     */
     sunperf_version_(&a, &b, &c);
     sprintf(b, "%d.%d.%d", a, b, c);
     return b;
#else
     return "unknown"
#fi
}

BEGIN_BENCH_DOC
BENCH_DOC("name", "sunperf")
BENCH_DOCF("version", mkvers)
END_BENCH_DOC


int can_do(struct problem *p)
{
     return (p->rank == 1 && problem_in_place(p));
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

void setup(struct problem *p)
{
     int n;
 
     BENCH_ASSERT(can_do(p));
     n = p->n[0];
 
     if (p->kind == PROBLEM_COMPLEX) {
	  WSAVE = bench_malloc((4 * n + 15) * sizeof(bench_real));
	  SINGLE_PRECISION ? cffti(n, WSAVE) : zffti(n, WSAVE);
     } else {
	  WSAVE = bench_malloc((2 * n + 15) * sizeof(bench_real));
	  SINGLE_PRECISION ? sffti(n, WSAVE) : dffti(n, WSAVE);
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
	       if (SINGLE_PRECISION)
		    for (i = 0; i < iter; ++i) {
			 cfftf(n, in, wsave);
		    }
	       else
		    for (i = 0; i < iter; ++i) {
			 zfftf(n, in, wsave);
		    }
	  } else {
	       if (SINGLE_PRECISION)
		    for (i = 0; i < iter; ++i) {
			 cfftb(n, in, wsave);
		    }
	       else
		    for (i = 0; i < iter; ++i) {
			 zfftb(n, in, wsave);
		    }
	  }
     } else {
#if 0
	  if (p->sign == -1) {
	       for (i = 0; i < iter; ++i) {
		    F77_FUNC(dfftf, DFFTF)(&n, in, wsave);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) {
		    F77_FUNC(dfftb, DFFTB)(&n, in, wsave);
	       }
	  }
#endif
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(WSAVE);
}
