/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "mpfun90")
BENCH_DOC("author", "David H. Bailey")
BENCH_DOC("year", "2001")
BENCH_DOC("language", "Fortran 90")
BENCH_DOC("email", "dhbailey@lbl.gov")
BENCH_DOC("url", "http://www.nersc.gov/~dhbailey/mpdist/mpdist.html")
BENCH_DOC("url-was-valid-on", "Mon Sep  2 00:48:38 EDT 2002")
BENCH_DOC("bibitem", 
	  "D. H. Bailey, J. of Supercomputing, p. 23-35 (March 1990)")
BENCH_DOC("notes", "Part of the MPFUN90 multi-precision arithmetic library.")
BENCH_DOC("notes", "Employs two-pass 4-step and Stockham FFT algorithms.")
BENCH_DOC("copyright",
"This software was written while the author was an employee of NASA.\n"
"This software has been approved by NASA for unrestricted distribution.\n"
"However, usage of this software is subject to the following:\n"
"\n"
"1. This software is offered without warranty of any kind, either expressed\n"
"   or implied.  The author would appreciate, however, any reports of bugs\n"
"   or other difficulties that may be encountered.\n"
"2. If modifications or enhancements to this software are made to this\n"
"   software by others, NASA Ames reserves the right to obtain this enhanced\n"
"   software at no cost and with no restrictions on its usage.\n"
"3. The author and NASA Ames are to be acknowledged in any published paper\n"
"   based on computations using this software.  Accounts of practical\n"
"   applications or other benefits resulting from this software are of\n"
"   particular interest.  Please send a copy of such papers to the author.\n")
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

/* We must call mpinix before calling any FFT routines, in order to
   initialize the trig tables. */
#define MPINIX_F77 F77_FUNC(mpinix, MPINIX)
extern void MPINIX_F77(unsigned int *n);

/* Bailey complex FFT: */
#define MPFFT1_F77 F77_FUNC(mpfft1, MPFFT1)
extern void MPFFT1_F77(int *is, int *m, int *n1, int *n2,
		       bench_complex *x, bench_complex *y);

/* Bailey real<->complex FFT routines.  These just wrap around the
   complex FFT with the usual n/2 trick, so I'm not sure they're worth
   benchmarking separately. */
#define MPFFTRC_F77 F77_FUNC(mpfftrc, MPFFTRC)
extern void MPFFTRC_F77(int *is, int *m, int *n, int *nsq,
			bench_complex *x, bench_complex *y);
#define MPFFTCR_F77 F77_FUNC(mpfftcr, MPFFTCR)
extern void MPFFTCR_F77(int *is, int *m, int *n, int *nsq,
			bench_complex *x, bench_complex *y);

static bench_complex *y = 0;
static int m, n1, n2;

static const int mpnsp1 = 2; /* must mirror definition in mpfun90.f!! */

void setup(struct problem *p)
{
     unsigned int n = p->n[0];
     BENCH_ASSERT(can_do(p));
     m = log_2(n);
     n1 = 1 << ((m + 1) / 2);
     n2 = 1 << (m - (m + 1) / 2);
     y = (bench_complex *) bench_malloc((n + n1 * mpnsp1) * 
					sizeof(bench_complex));
     MPINIX_F77(&n);
}

void doit(int iter, struct problem *p)
{
     int is = (p->sign > 0) ? 1 : -1;
     bench_complex *in = (bench_complex *) p->in;
     int i;

     for (i = 0; i < iter; ++i) {
	  MPFFT1_F77(&is, &m, &n1, &n2, in, y);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(y);
}
