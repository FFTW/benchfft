/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "fftpack")
BENCH_DOC("author", "Paul N. Swarztrauber")
BENCH_DOC("year", "1985")
BENCH_DOC("version", "4")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("url", "ftp://ftp.netlib.org:/fftpack")
BENCH_DOC("url-was-valid-on", "Thu Jul 12 20:26:24 EDT 2001")
BENCH_DOC("bibitem", 
	  "P. N. Swarztrauber, Vectorizing the FFTs, Parallel Computations, p. 51-83 (1982).")
BENCH_DOC("copyright", "Public domain")
END_BENCH_DOC


int can_do(struct problem *p)
{
     return (SINGLE_PRECISION && p->rank == 1 && problem_in_place(p));
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

extern void F77_FUNC(cffti, CFFTI)();
extern void F77_FUNC(rffti, RFFTI)();
extern void F77_FUNC(cfftf, CFFTF)();
extern void F77_FUNC(rfftf, RFFTF)();
extern void F77_FUNC(cfftb, CFFTB)();
extern void F77_FUNC(rfftb, RFFTB)();

void setup(struct problem *p)
{
     int n;
     const int extra_locations = 500;  
     /* not clear how many locations are actually needed.  FFTPACK
	plays it dirty by confusing 2int = 1 double, which is false
	on 64bit machines */
 
     BENCH_ASSERT(can_do(p));
     n = p->n[0];
 
     if (p->kind == PROBLEM_COMPLEX) {
	  WSAVE = bench_malloc((4 * n + 15 + extra_locations)
			       * sizeof(bench_real));
	  F77_FUNC(cffti, CFFTI)(&n, WSAVE);
     } else {
	  WSAVE = bench_malloc((2 * n + 15 + extra_locations)
			       * sizeof(bench_real));
	  F77_FUNC(rffti, RFFTI)(&n, WSAVE);
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
		    F77_FUNC(cfftf, CFFTF)(&n, in, wsave);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) {
		    F77_FUNC(cfftb, CFFTB)(&n, in, wsave);
	       }
	  }
     } else {
	  if (p->sign == -1) {
	       for (i = 0; i < iter; ++i) {
		    F77_FUNC(rfftf, RFFTF)(&n, in, wsave);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) {
		    F77_FUNC(rfftb, RFFTB)(&n, in, wsave);
	       }
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(WSAVE);
}
