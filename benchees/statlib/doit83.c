/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "as83")
BENCH_DOC("package", "StatLib")
BENCH_DOC("author", "Donald M. Monro")
BENCH_DOC("year", "1975")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("url", "http://lib.stat.cmu.edu/apstat/83")
BENCH_DOC("url-was-valid-on", "Sun Jan 13 14:14:59 EST 2002")
BENCH_DOC("notes", "Downloaded from the StatLib repository at Carnegie Mellon University, in the Applied Statistics algorithm collection.")
BENCH_DOC("notes", "Complex data are stored in separate real/imag arrays.")
BENCH_DOC("notes", "The forward transform is scaled.")
BENCH_DOC("notes", "The inverse transform is computed via the inverse transform by a pair of complex conjugations.")
BENCH_DOC("bibitem", 
	  "Donald M. Monro, Algorithm AS 83: Complex discrete fast Fourier transform, Appl. Statist. 24 (1), pp. 153-160 (1975).")
BENCH_DOC("copyright", "The Royal Statistical Society holds the copyright to these routines, but has given its permission for their distribution provided that no fee is charged.")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (SINGLE_PRECISION &&
	     p->rank == 1 &&
	     p->n[0] >= 4 && p->n[0] <= (1<<20) &&
	     problem_complex_power_of_two(p, 1)
	  );
}

#define FASTF F77_FUNC(fastf,FASTF)
extern void FASTF(bench_real *re, bench_real *im,
		  int *n, int *isn);

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
}

void doit(int iter, struct problem *p)
{
     bench_real *in_i, *in_r = (bench_real *) p->in;
     int i;
     int n = p->n[0];
     int isn = (p->sign > 0) ? -1 : 1;
     in_i = in_r + n;

     for (i = 0; i < iter; ++i) {
	  FASTF(in_r, in_i, &n, &isn);
     }
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, -1);
}

void copy_c2c_from(struct problem *p, bench_complex *in)
{
     bench_real *x = (bench_real *) p->in;
     copy_c2ri(in, x, x + p->size, p->size);
}

void copy_c2c_to(struct problem *p, bench_complex *out)
{
     bench_real *x = (bench_real *) p->out;
     copy_ri2c(x, x + p->size, out, p->size);
}

void done(struct problem *p)
{
     UNUSED(p);
}
