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
END_BENCH_DOC

#ifdef BENCHFFT_SINGLE
#define FFT2I cfft2i
#define FFT2F cfft2f
#define FFT2B cfft2b
#define RFFT2I rfft2i
#define RFFT2F rfft2f
#define RFFT2B rfft2b
#else
#define FFT2I zfft2i
#define FFT2F zfft2f
#define FFT2B zfft2b
#define RFFT2I dfft2i
#define RFFT2F dfft2f
#define RFFT2B dfft2b
#endif

int can_do(struct problem *p)
{
     return (p->rank == 2 && 
	     /* real problems can be out of place */
	     (p->kind == PROBLEM_REAL || problem_in_place(p)));
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_unpacked(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_unpacked(p, in, -1.0);
}


void copy_r2c(struct problem *p, bench_complex *out)
{
     copy_r2c_unpacked(p, out);	  
}

void copy_c2r(struct problem *p, bench_complex *in)
{
     copy_c2r_unpacked(p, in);
}

static void *WSAVE;
static int LWORK;

void setup(struct problem *p)
{
     int n, m;
 
     BENCH_ASSERT(can_do(p));
     m = p->n[1];
     n = p->n[0];
 
     if (p->kind == PROBLEM_COMPLEX) {
	  LWORK = 4 * (n + m) + 30; /* the value 2 * (m + n) + 30 in the
				       manual appears to be bogus */
	  WSAVE = bench_malloc(LWORK * sizeof(bench_real));
	  FFT2I(m, n, WSAVE);
     } else {
	  LWORK = 2 * m + 3 * n + 30;
	  WSAVE = bench_malloc(LWORK * sizeof(bench_real));
	  RFFT2I(m, n, WSAVE);
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     int m = p->n[1];
     int n = p->n[0];
     void *in = p->in;
     void *wsave = WSAVE;
     int lwork = LWORK;

     if (p->kind == PROBLEM_COMPLEX) {
	  int lda = m;
	  if (p->sign == -1) {
	    for (i = 0; i < iter; ++i) 
	      FFT2F(m, n, in, lda, wsave, lwork);
	  } else {
	    for (i = 0; i < iter; ++i) 
	      FFT2B(m, n, in, lda, wsave, lwork);
	  }
     } else {
	  if (p->sign == -1) {
	       void *out = p->out;
	       char place = p->in_place ? 'i' : 'o';
	       int lda = 2 * (1 + m / 2);
	       int ldb = lda / 2;
	       for (i = 0; i < iter; ++i) {
		    RFFT2F(place, 0, m, n, in, lda, out, ldb, wsave, lwork);
	       }
	  } else {
	       int lda = 2 * (1 + m / 2);
	       int ldb;
	       void *out;
	       char place;

	       if (p->out == in) {
		 place = 'i'; ldb = lda; out = in;
	       } else {
		 place = 'o'; out = in; in = p->out; ldb = lda / 2;
	       }
	       for (i = 0; i < iter; ++i) {
		    RFFT2B(place, m, n, in, lda, out, ldb, wsave, lwork);
 	       }
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(WSAVE);
}
