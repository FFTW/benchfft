/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "mpfun77")
BENCH_DOC("author", "David H. Bailey")
BENCH_DOC("year", "1988")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("email", "dhbailey@lbl.gov")
BENCH_DOC("url", "http://www.nersc.gov/~dhbailey/mpdist/mpdist.html")
BENCH_DOC("url-was-valid-on", "Mon Sep  2 00:48:38 EDT 2002")
BENCH_DOC("bibitem", 
  "D. H. Bailey, Intl. J. of Supercomp. Appl., p. 82-87 (Spring 1988).")
BENCH_DOC("notes", "Complex data are stored in separate real/imag arrays.")
BENCH_DOC("notes", "Part of the MPFUN multi-precision arithmetic library.")
BENCH_DOC("notes", "Employs 4-step and Stockham FFT algorithms.")
BENCH_DOC("notes", "Patched by S. G. Johnson to allow arrays to be allocated from C caller.")
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
   initialize the trig tables.  Getting this to work, however,
   requires slight hacking to separate mpinix's "u" array into an
   external parameter that we can allocate. */

#define MPINIX_F77 F77_FUNC(mpinix, MPINIX)
extern void MPINIX_F77(int *m, bench_real *u);

/* Bailey complex FFT: */
#define MPCFFT_F77 F77_FUNC(mpcfft, MPCFFT)
extern void MPCFFT_F77(int *is, int *m, bench_real *x, bench_real *y,
		       bench_real *u);

/* Bailey real<->complex FFT routines.  These just wrap around the
   complex FFT with the usual n/2 trick, so I'm not sure they're worth
   benchmarking separately. */
#define MPRCFT_F77 F77_FUNC(mprcft, MPRCFT)
extern void MPRCFT_F77(int *is, int *m, bench_real *x, bench_real *y,
		       bench_real *u);
#define MPCRFT_F77 F77_FUNC(mpcrft, MPCRFT)
extern void MPCRFT_F77(int *is, int *m, bench_real *x, bench_real *y,
		       bench_real *u);

static bench_real *y = 0, *u = 0;
static int m;

void setup(struct problem *p)
{
     unsigned int n = p->n[0];
     BENCH_ASSERT(can_do(p));
     y = (bench_real *) bench_malloc(n * 2 * sizeof(bench_real));
     u = (bench_real *) bench_malloc(n * 8 * sizeof(bench_real));
     m = log_2(n);
     MPINIX_F77(&m, u);
}

void problem_ccopy_from(struct problem *p, bench_complex *in)
{
     bench_real *x = (bench_real *) p->in;
     copy_c2ri(in, x, x + p->n[0], p->n[0]);
}

void problem_ccopy_to(struct problem *p, bench_complex *out)
{
     bench_real *x = (bench_real *) p->in;
     copy_ri2c(x, x + p->n[0], out, p->n[0]);
}

void doit(int iter, struct problem *p)
{
     int is = (p->sign > 0) ? 1 : -1;
     bench_real *in = (bench_real *) p->in;
     int i;

     for (i = 0; i < iter; ++i) {
	  MPCFFT_F77(&is, &m, in, y, u);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(u);
     bench_free(y);
}
