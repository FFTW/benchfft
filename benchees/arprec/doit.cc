/* this program is in the public domain */

#include "bench-user.h"
#include "mp/mpreal.h"

BEGIN_BENCH_DOC
BENCH_DOC("name", "arprec")
BENCH_DOC("author", "David H. Bailey")
BENCH_DOC("email", "dhbailey@lbl.gov")
BENCH_DOC("author", "Yozo Hida")
BENCH_DOC("author", "Xiaoye S. Li")
BENCH_DOC("author", "Brandon Thompson")
BENCH_DOC("year", "2002")
BENCH_DOC("version", "2002-08-27")
BENCH_DOC("language", "C++")
BENCH_DOC("url", "http://crd.lbl.gov/~dhbailey/mpdist/")
BENCH_DOC("url-was-valid-on", "Thu Aug 21 17:07:36 EDT 2003")
BENCH_DOC("bibitem", 
	  "D. H. Bailey, J. of Supercomputing, p. 23-35 (March 1990)")
BENCH_DOC("notes", "Part of the ARPREC arbitrary-precision arithmetic library.")
BENCH_DOC("notes", "Employs two-pass 4-step and Stockham FFT algorithms.")
BENCH_DOC("notes", "Based on code in DHB's MPFUN package")
BENCH_DOC("copyright",
	  "Copyright (c) 2002\n"
	  "Revised  2002-08-27\n"
	  "\n"
	  "This software is provided for research purposes only.\n"
	  "Commercial usage requires license agreement.\n")
BENCH_DOC("notes", "A new version was released on 2003-07-14, but has not made it into our benchmark yet.  Sorry about that.")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (DOUBLE_PRECISION &&
	     p->rank == 1 &&
	     p->kind == PROBLEM_COMPLEX &&
	     p->n[0] >= 4 &&
	     problem_power_of_two(p, 1)
	  );
}

/* There are also real<->complex FFT routines.  These just wrap around the
   complex FFT with the usual n/2 trick, so I'm not sure they're worth
   benchmarking separately. */

static double *y = 0;
static int m, m1, m2, n1, n2;

void setup(struct problem *p)
{
     unsigned int n = p->n[0];
     BENCH_ASSERT(can_do(p));
     m = log_2(n);
     n1 = 1 << (m1 = ((m + 1) / 2));
     n2 = 1 << (m2 = (m - (m + 1) / 2));
     y = (double *) bench_malloc((n + n1 * mp::mpnsp1) * 
					sizeof(bench_complex));
     mp_real::mpinix(n);
}

void doit(int iter, struct problem *p)
{
     int is = p->sign;
     double *in = (double *) p->in;
     int i;

     for (i = 0; i < iter; ++i) {
	  mp_real::mpfft1(is, m, m1, m2, in, y);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(y);
}
