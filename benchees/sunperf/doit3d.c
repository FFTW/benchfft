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
     return "unknown"
#endif
}

BEGIN_BENCH_DOC
BENCH_DOC("name", "sunperf")
BENCH_DOCF("version", mkvers)
BENCH_DOC("notes", "The library has 3D real transforms but I can't get them to work")
END_BENCH_DOC

#ifdef BENCHFFT_SINGLE
#define FFT3I cfft3i
#define FFT3F cfft3f
#define FFT3B cfft3b
#define RFFT3I rfft3i
#define RFFT3F rfft3f
#define RFFT3B rfft3b
#else
#define FFT3I zfft3i
#define FFT3F zfft3f
#define FFT3B zfft3b
#define RFFT3I dfft3i
#define RFFT3F dfft3f
#define RFFT3B dfft3b
#endif

int can_do(struct problem *p)
{
     return (p->rank == 3 
	     && p->kind == PROBLEM_COMPLEX /* can't get real to work properly */
	     /* real problems can be out of place */
	     && (p->kind == PROBLEM_REAL || problem_in_place(p)));
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
     int n, m, k;
 
     BENCH_ASSERT(can_do(p));
     m = p->n[2];
     n = p->n[1];
     k = p->n[0];
 
     if (p->kind == PROBLEM_COMPLEX) {
	  LWORK = 4 * (n + m + k) + 45; /* the value 2 * (m + n) + 45 in the
					   manual appears to be bogus */
	  WSAVE = bench_malloc(LWORK * sizeof(bench_real));
	  FFT3I(m, n, k, WSAVE);
     } else {
          /* it's easier to win the lottery than to get this number
	     right.  Each manpage shows a different number */
	  LWORK = 4 * (m + n + k) + 45;
	  WSAVE = bench_malloc(LWORK * sizeof(bench_real));
	  RFFT3I(m, n, k, WSAVE);
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     int m = p->n[2];
     int n = p->n[1];
     int k = p->n[0];
     void *in = p->in;
     void *wsave = WSAVE;
     int lwork = LWORK;

     if (p->kind == PROBLEM_COMPLEX) {
	  int lda = m;
	  int ld2a = n;
	  if (p->sign == -1) {
	    for (i = 0; i < iter; ++i) 
	      FFT3F(m, n, k, in, lda, ld2a, wsave, lwork);
	  } else {
	    for (i = 0; i < iter; ++i) 
	      FFT3B(m, n, k, in, lda, ld2a, wsave, lwork);
	  }
     } else {
	  if (p->sign == -1) {
	       void *out = p->out;
 	       char place = in == out ? 'i' : 'o';
	       int lda = 2 * (1 + m / 2);
	       int ldb = lda / 2;
	       for (i = 0; i < iter; ++i) {
		    RFFT3F(place, 0, m, n, k, in, lda, out, ldb, wsave, lwork);
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
		    RFFT3B(place, m, n, k, in, lda, out, ldb, wsave, lwork);
 	       }
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(WSAVE);
}
