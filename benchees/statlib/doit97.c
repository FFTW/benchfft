/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "as97")
BENCH_DOC("package", "StatLib")
BENCH_DOC("author", "Donald M. Monro")
BENCH_DOC("year", "1976")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("url", "http://lib.stat.cmu.edu/apstat/97")
BENCH_DOC("url-was-valid-on", "Sun Jan 13 14:14:59 EST 2002")
BENCH_DOC("notes", "Downloaded from the StatLib repository at Carnegie Mellon University, in the Applied Statistics algorithm collection.")
BENCH_DOC("notes", "The forward transform is scaled.")
BENCH_DOC("bibitem", 
	  "Donald M. Monro, Algorithm AS 97: Real discrete fast Fourier transform, Appl. Statist. 25 (2), pp. 166-172 (1976).")
BENCH_DOC("copyright", "The Royal Statistical Society holds the copyright to these routines, but has given its permission for their distribution provided that no fee is charged.")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (SINGLE_PRECISION &&
	     p->rank == 1 &&
	     p->n[0] >= 8 && p->n[0] <= (1<<21) &&
	     problem_real_power_of_two(p, 1)
	  );
}

#define FORRT F77_FUNC(forrt,FORRT)
extern void FORRT(bench_real *x, int *n);

#define REVRT F77_FUNC(revrt,REVRT)
extern void REVRT(bench_real *x, int *n);

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
}

void doit(int iter, struct problem *p)
{
     bench_real *x = (bench_real *) p->in;
     int n = p->n[0];
     int i;
     
     if (p->sign < 0) {
	  for (i = 0; i < iter; ++i) {
	       FORRT(x, &n);
	  }
     }
     else {
	  for (i = 0; i < iter; ++i) {
	       REVRT(x, &n);
	  }
     }
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, -1);
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     int i, n = p->size;
     bench_real *x = (bench_real *) p->out;
     
     for (i = n/2 + 1; i < n - (i - n/2); ++i) {
	  bench_real y = x[i];
	  x[i] = x[n - (i - n/2)];
	  x[n - (i - n/2)] = y;
     }
     
     copy_h2c_1d_halfcomplex(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     int i, n = p->size;
     bench_real *x = (bench_real *) p->in;
     
     copy_c2h_1d_halfcomplex(p, in, -1.0);

     for (i = n/2 + 1; i < n - (i - n/2); ++i) {
	  bench_real y = x[i];
	  x[i] = x[n - (i - n/2)];
	  x[n - (i - n/2)] = y;
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
