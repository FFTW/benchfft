/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "temperton")
BENCH_DOC("author", "Clive Temperton, modified by Russ Rew")
BENCH_DOC("year", "1980")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("bibitem", 
	  "The package was written by Clive Temperton at ECMWF in\n"
	  "November, 1978.  It was modified, documented, and tested\n"
	  "for NCAR by Russ Rew in September, 1980.\n")
BENCH_DOC("url", "http://www.caps.ou.edu/ARPS/")
BENCH_DOC("url-was-valid-on", "Mon Jul 16 00:03:31 EDT 2001")
BENCH_DOC("notes", "Can be found online in the ARPS weather-simulator")
BENCH_DOC("notes", "The forward (real->complex) transform is scaled.")
END_BENCH_DOC


int can_do(struct problem *p)
{
     return (SINGLE_PRECISION &&
	     p->rank == 1 &&
	     p->kind == PROBLEM_REAL &&
	     check_prime_factors(p->n[0], 5) &&
	     p->n[0] % 2 == 0 &&
	     p->n[0] > 4 &&
	     problem_in_place(p)
	     /* && if n has more than one square-free factor, the
		   product of the square-free factors must be .le. 210 
		&& n has at most 15 prime factors
	     */
	  );
}

#define FFT99_F77 F77_FUNC(fft99,FFT99)
#define FFT991_F77 F77_FUNC(fft991,FFT991)
#define SET99_F77 F77_FUNC(set99,SET99)

extern void FFT99_F77(bench_real *a, bench_real *work,
		      bench_real *trigs, int ifax[13],
		      int *inc, int *jump, int *n, int *lot, int *isign);
extern void FFT991_F77(bench_real *a, bench_real *work,
		       bench_real *trigs, int ifax[13],
		       int *inc, int *jump, int *n, int *lot, int *isign);

extern void SET99_F77(bench_real *trigs, int ifax[13], int *n);

bench_real *work;
bench_real *trigs;
int ifax[13];

void setup(struct problem *p)
{
     int n = p->n[0];

     BENCH_ASSERT(can_do(p));
     trigs = (bench_real *) bench_malloc((3*n/2 + 1) * sizeof(bench_real));
     SET99_F77(trigs, ifax, &n);
     work = (bench_real *) bench_malloc((n + 1) * sizeof(bench_real));
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_unpacked(p, out, -1);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_unpacked(p, in, -1.0);
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, -1.0);
}

void after_problem_ccopy_from(struct problem *p, bench_complex *in)
{
     if (p->sign == -1) {
	  int n = p->n[0];
	  bench_real *a = (bench_real *) p->in;
	  a[n] = a[n+1] = 0.0;
     }
}

void doit(int iter, struct problem *p)
{
     int i, n = p->n[0];
     int inc = 1, jump = 1, lot = 1, isign = p->sign;
     bench_real *a = (bench_real *) p->in;

     for (i = 0; i < iter; ++i) {
	  FFT991_F77(a, work, trigs, ifax, &inc, &jump, &n, &lot, &isign);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(work);
     bench_free(trigs);
}
