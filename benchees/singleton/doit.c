/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "singleton")
BENCH_DOC("author", "R. C. Singleton")
BENCH_DOC("year", "1968")
BENCH_DOC("language", "FORTRAN")
END_BENCH_DOC


int can_do(struct problem *p)
{
     return (
	  sizeof(bench_real) == sizeof(float) && 
	  p->rank == 1 &&
	  p->kind == PROBLEM_COMPLEX &&
	  p->in == p->out &&
	  check_prime_factors(p->n[0], 23) &&
	  problem_in_place(p)
	  /* &&  if n has more than one square-free factor, the product of the
	         square-free factors must be .le. 210
	     &&  n has at most 15 prime factors
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
     int n = p->n[0];
     int ntot = n;
     int nspan = n;
     int isn = (p->sign > 0) ? 2 : -2;
     bench_complex *in = p->in;
     int i;

     for (i = 0; i < iter; ++i) {
	  F77_FUNC(fft, FFT)(&c_re(in[0]), &c_im(in[0]),
			     &ntot, &n, &nspan, &isn);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
