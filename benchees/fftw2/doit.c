/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

#if defined(BENCHFFT_SINGLE) && defined(HAVE_SFFTW_H) && defined(HAVE_SRFFTW_H)
#  include <sfftw.h>
#  include <srfftw.h>
#else
#  include <fftw.h>
#  include <rfftw.h>
#endif

static const char *mkvers(void)
{
     return fftw_version;
}

BEGIN_BENCH_DOC
BENCH_DOC("name", NAME)
BENCH_DOCF("version", mkvers)
BENCH_DOC("year", "2003")
BENCH_DOC("author", "Matteo Frigo")
BENCH_DOC("author", "Steven G. Johnson")
BENCH_DOC("email", "fftw@fftw.org")
BENCH_DOC("url", "http://www.fftw.org")
BENCH_DOC("url-was-valid-on", "Fri Mar 28 18:46:22 EST 2003")
BENCH_DOC("language", "C")
BENCH_DOC("language", "Objective Caml")
BENCH_DOC("copyright",
"Copyright (c) 1997-1999, 2003 Massachusetts Institute of Technology\n"
"\n"
"This program is free software; you can redistribute it and/or modify\n"
"it under the terms of the GNU General Public License as published by\n"
"the Free Software Foundation; either version 2 of the License, or\n"
"(at your option) any later version.\n"
"\n"
"This program is distributed in the hope that it will be useful,\n"
"but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"GNU General Public License for more details.\n"
"\n"
"You should have received a copy of the GNU General Public License\n"
"along with this program; if not, write to the Free Software\n"
"Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA\n")
BENCH_DOC("bibitem", "M. Frigo and S. G. Johnson, FFTW: An adaptive software architecture for the FFT, Proc. ICASSP 3, 1381-1384 (1998)")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (sizeof(fftw_real) == sizeof(bench_real) && (p->rank == 1));
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_1d_halfcomplex(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_halfcomplex(p, in, -1.0);
}

static void *plan;

void setup(struct problem *p)
{
     int n;
 
     BENCH_ASSERT(can_do(p));
     n = p->n[0];
 
     if (p->kind == PROBLEM_COMPLEX) {
	  plan = fftw_create_plan_specific(
	       n, p->sign,
	       (PLAN_FLAGS | (problem_in_place(p) ? FFTW_IN_PLACE : 0)),
	       p->in, 1,
	       p->out, 1);
     } else {
	  plan = rfftw_create_plan_specific(
	       n, p->sign,
	       (PLAN_FLAGS | (problem_in_place(p) ? FFTW_IN_PLACE : 0)),
	       p->in, 1,
	       p->out, 1);
     }
     BENCH_ASSERT(plan);
}

void doit(int iter, struct problem *p)
{
     int i;
     void *in = p->in;
     void *out = p->out;
     void *PLAN = plan;

     if (out == in)
	  out = 0; /* scratch space for in-place */

     if (p->kind == PROBLEM_COMPLEX) {
	  for (i = 0; i < iter; ++i) {
	       fftw_one(PLAN, in, out);
	  }
     } else {
	  for (i = 0; i < iter; ++i) {
	       rfftw_one(PLAN, in, out);
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     if (p->kind == PROBLEM_COMPLEX) {
	  fftw_destroy_plan(plan);
     } else {
	  rfftw_destroy_plan(plan);
     }
}
