/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "singleton")
BENCH_DOC("author", "R. C. Singleton")
BENCH_DOC("year", "1968")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("url", "http://www.netlib.org/go/fft.f")
BENCH_DOC("url", "http://www.netlib.org/go/realtr.f")
BENCH_DOC("url-was-valid-on", "Thu Jul 12 20:26:24 EDT 2001")
BENCH_DOC("notes", "Used complex data format (alternating real/imag. parts), which usually seems to be faster than separate real/imag. arrays.")
BENCH_DOC("bibitem", 
	  "R. C. Singleton, An algorithm for computing the mixed radix fast Fourier transform, IEEE Trans. on Audio and Electroacoustics AU-17, no. 2, p. 93-103 (June, 1969).")
END_BENCH_DOC


int can_do(struct problem *p)
{
     return (SINGLE_PRECISION &&
	     p->rank == 1 &&
	     (p->kind == PROBLEM_COMPLEX || 
	      (p->kind == PROBLEM_REAL && p->n[0] % 2 == 0)) &&
	     check_prime_factors(p->n[0], 23) &&
	     problem_in_place(p)
	     /* && if n has more than one square-free factor, the
		   product of the square-free factors must be .le. 210 
		&& n has at most 15 prime factors
	     */
	  );
}

#define FFT F77_FUNC(fft,FFT)
extern void FFT(bench_real *re, bench_real *im,
		int *ntot, int *n, int *nspan, int *isn);

#define REALTR F77_FUNC(realtr,REALTR)
extern void REALTR(bench_real *re, bench_real *im,
		   int *n, int *isn);

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
}

void doit(int iter, struct problem *p)
{
     bench_complex *in = p->in;
     int i;

     if (p->kind == PROBLEM_COMPLEX) {
	  int n = p->n[0];
	  int ntot = n;
	  int nspan = n;
	  int isn = (p->sign > 0) ? 2 : -2;

	  for (i = 0; i < iter; ++i) {
	       FFT(&c_re(in[0]), &c_im(in[0]), &ntot, &n, &nspan, &isn);
	  }
     }
     else {
	  int n = p->n[0] / 2;
	  int ntot = n;
	  int nspan = n;
	  int isn = (p->sign > 0) ? -2 : 2;

	  if (p->sign < 0)
	       for (i = 0; i < iter; ++i) {
		    FFT(&c_re(in[0]), &c_im(in[0]), &ntot, &n, &nspan, &isn);
		    REALTR(&c_re(in[0]), &c_im(in[0]), &n, &isn);
	       }
	  else
	       for (i = 0; i < iter; ++i) {
		    REALTR(&c_re(in[0]), &c_im(in[0]), &n, &isn);
		    FFT(&c_re(in[0]), &c_im(in[0]), &ntot, &n, &nspan, &isn);
	       }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_unpacked(p, out, 1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_unpacked(p, in, 1.0);
}

/* Forwards real transform is scaled by 2.0 */
void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     if (p->kind == PROBLEM_REAL && p->sign < 0) {
	  bench_complex x;
	  c_re(x) = 0.5;
	  c_im(x) = 0;
	  
	  cascale(out, p->size, x);
     }
}
