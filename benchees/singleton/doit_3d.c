/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "singleton-3d")
BENCH_DOC("author", "R. C. Singleton")
BENCH_DOC("year", "1968")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("url", "http://www.netlib.org/go/fft.f")
BENCH_DOC("url-was-valid-on", "Thu Jul 12 20:26:24 EDT 2001")
BENCH_DOC("bibitem", 
	  "R. C. Singleton, An algorithm for computing the mixed radix fast Fourier transform, IEEE Trans. on Audio and Electroacoustics AU-17, no. 2, p. 93-103 (June, 1969).")
END_BENCH_DOC


int can_do(struct problem *p)
{
     return (SINGLE_PRECISION &&
	     p->rank == 3 &&
	     p->kind == PROBLEM_COMPLEX &&
	     problem_in_place(p) &&
	     check_prime_factors(p->n[0], 23) &&
	     check_prime_factors(p->n[1], 23) &&
	     check_prime_factors(p->n[2], 23) &&
	     problem_in_place(p)
	     /* && if n has more than one square-free factor, the
		   product of the square-free factors must be .le. 210 
		&& n has at most 15 prime factors
	     */
	  );
}

extern void F77_FUNC(fft, FFT)();

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
}

void doit(int iter, struct problem *p)
{
     int n3 = p->n[0];
     int n2 = p->n[1];
     int n1 = p->n[2];
     int ntot = n1*n2*n3;
     int nspan12 = n1*n2;
     int isn = (p->sign > 0) ? 2 : -2;
     bench_complex *in = p->in;
     int i;

     for (i = 0; i < iter; ++i) {
	  F77_FUNC(fft, FFT)(&c_re(in[0]), &c_im(in[0]),
			     &ntot, &n1, &n1, &isn);
	  F77_FUNC(fft, FFT)(&c_re(in[0]), &c_im(in[0]),
			     &ntot, &n2, &nspan12, &isn);
	  F77_FUNC(fft, FFT)(&c_re(in[0]), &c_im(in[0]),
			     &ntot, &n3, &ntot, &isn);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
