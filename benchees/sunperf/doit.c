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
     static char buf[128];

     /* 
	using FORTRAN interface.  The documented C interface cannot
	possibly work because it takes a, b, c by value and it returns
	void
     */
     sunperf_version_(&a, &b, &c);
     sprintf(buf, "%d.%d.%d", a, b, c);
     return buf;
#else
     return "unknown";
#endif
}

BEGIN_BENCH_DOC
BENCH_DOC("name", "sunperf")
BENCH_DOCF("version", mkvers)
BENCH_DOC("package", "Sun Performance Library (SUNPERF)")
BENCH_DOC("notes", "1d transforms")
END_BENCH_DOC

#ifdef BENCHFFT_SINGLE
#define FFTI cffti
#define FFTF cfftf
#define FFTB cfftb
#define RFFTI rffti
#define RFFTF rfftf
#define RFFTB rfftb
#else
#define FFTI zffti
#define FFTF zfftf
#define FFTB zfftb
#define RFFTI dffti
#define RFFTF dfftf
#define RFFTB dfftb
#endif

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

static void *WSAVE;

void setup(struct problem *p)
{
     int n;
 
     BENCH_ASSERT(can_do(p));
     n = p->n[0];
 
     if (p->kind == PROBLEM_COMPLEX) {
	  WSAVE = bench_malloc((4 * n + 15) * sizeof(bench_real));
	  FFTI(n, WSAVE);
     } else {
	  WSAVE = bench_malloc((2 * n + 15) * sizeof(bench_real));
	  RFFTI(n, WSAVE);
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     int n = p->n[0];
     void *in = p->in;
     void *wsave = WSAVE;

     if (p->kind == PROBLEM_COMPLEX) {
	  if (p->sign == -1) {
	    for (i = 0; i < iter; ++i) 
	      FFTF(n, in, wsave);
	  } else {
	    for (i = 0; i < iter; ++i) 
	      FFTB(n, in, wsave);
	  }
     } else {
	  if (p->sign == -1) {
	       for (i = 0; i < iter; ++i) {
		    RFFTF(n, in, wsave);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) {
		    RFFTB(n, in, wsave);
	       }
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(WSAVE);
}
