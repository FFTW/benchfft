/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "napack")
BENCH_DOC("package", "NAPACK")
BENCH_DOC("author", "William W. Hager")
BENCH_DOC("year", "1987")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("url", "http://www.netlib.org/napack/")
BENCH_DOC("url-was-valid-on", "Fri Aug 17 23:36:01 EDT 2001")
BENCH_DOC("bibitem", 
	  "William W. Hager, Applied Numerical Linear Algebra (Prentice-Hall, 1987).")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (SINGLE_PRECISION &&
	     p->rank == 1 &&
	     p->kind == PROBLEM_COMPLEX &&
	     /* check_prime_factors(p->n[0], 23) && */
             p->n[0] > 1 && /* wrong for n=1 */
	     problem_in_place(p)
	  );
}

#define FFT F77_FUNC(fft,FFT)
extern void FFT(bench_complex *a, int *n, bench_complex *work);

#define FFC F77_FUNC(ffc,FFC)
extern void FFC(bench_complex *a, int *n, bench_complex *work);

static bench_complex *work = 0;

void setup(struct problem *p)
{
     unsigned int i;
     BENCH_ASSERT(can_do(p));
     work = (bench_complex *) bench_malloc(sizeof(bench_complex) * p->size);
}

void doit(int iter, struct problem *p)
{
     bench_complex *in = p->in;
     int n = p->size;
     int i;

     if (p->sign > 0)
	  for (i = 0; i < iter; ++i) {
	       FFT(in, &n, work);
	  }
     else
	  for (i = 0; i < iter; ++i) {
	       FFC(in, &n, work);
	  }
}

void done(struct problem *p)
{
     bench_free(work);
}
