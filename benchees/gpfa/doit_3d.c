/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "gpfa-3d")
BENCH_DOC("author", "Clive Temperton")
BENCH_DOC("year", "1992 (?)")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("bibitem", 
	  "C. Temperton, A Generalized Prime Factor Fft Algorithm "
	  "For Any N = (2**P)(3**Q)(5**R), SIAM J. Sci. Stat. Comp. 13 (3),"
	  "p. 676-686 (May 1992).")
END_BENCH_DOC

static const unsigned int NMAX = 256; /* must match constant in gpfft3.f */

int can_do(struct problem *p)
{
     return (SINGLE_PRECISION && 
	     p->rank == 3 &&
	     p->kind == PROBLEM_COMPLEX &&
	     p->n[0] <= NMAX/2 &&
	     p->n[1] <= NMAX/2 &&
	     p->n[2] <= NMAX/2 &&
	     problem_in_place(p) &&
	     check_prime_factors(p->n[0], 5) &&
	     check_prime_factors(p->n[1], 5) &&
	     check_prime_factors(p->n[2], 5));
}

extern void F77_FUNC(gpf3d, GPF3D)(bench_complex *c, unsigned int *id,
				   unsigned int *nn, int *is);

unsigned int nn[3];

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
     nn[0] = p->n[2];
     nn[1] = p->n[1];
     nn[2] = p->n[0];
     doit(1, p); /* FFT initializes trig tables on first call */
}

void doit(int iter, struct problem *p)
{
     int i;
     bench_complex *in = p->in;
     int isign = -p->sign;

     for (i = 0; i < iter; ++i) {
	  F77_FUNC(gpf3d, GPF3D)(in, nn, nn, &isign);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
