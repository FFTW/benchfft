/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

#define FFTCACHE 16
#define STRINGIFYx(X) #X
#define STRINGIFY(X) STRINGIFYx(X)

BEGIN_BENCH_DOC
BENCH_DOC("name", "goedecker")
BENCH_DOC("author", "Stefan Goedecker")
BENCH_DOC("year", "1993")
BENCH_DOC("url", "http://www.abinit.org/")
BENCH_DOC("url-was-valid-on", "Thu May 22 00:37:07 EDT 2003")
BENCH_DOC("copyright",
	  "Copyright by Stefan Goedecker, Ithaca, NY USA, July 14, 1993\n"
	  "Copyright (C) 1998-2001 ABINIT group (DCA, XG)\n"
	  "This file is distributed under the terms of the\n"
	  "GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .")
BENCH_DOC("notes",
	  "This routine was downloaded as part of the ABINIT-3.1.3 software.\n")
BENCH_DOC("notes",
	  "Slightly modified for the benchmark.  See sgfft.tar.gz for the original code.")
BENCH_DOC("notes",
	  "The fftcache parameter is set to " STRINGIFY(FFTCACHE) "KB.  This may be wrong on your system.")
END_BENCH_DOC

#define SG_FFT_ F77_FUNC_(sg_fft, SG_FFT)
extern void SG_FFT_();
 
int can_do(struct problem *p)
{
     return (p->rank == 3 &&
	     DOUBLE_PRECISION &&
	     p->kind == PROBLEM_COMPLEX &&
	     p->n[0] <= 2048 && p->n[1] <= 2048 && p->n[2] <= 2048 &&
	     check_prime_factors(p->size, 5) &&
	     !problem_in_place(p));
}

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
     /* nothing to do */
}

void doit(int iter, struct problem *p)
{
     int i;
     int n1 = p->n[2];
     int n2 = p->n[1];
     int n3 = p->n[0];
     int fftcache = FFTCACHE;
     void *in = p->in;
     void *out = p->out;
     double sign = p->sign;

     for (i = 0; i < iter; ++i) {
	  SG_FFT_(&fftcache, &n1, &n2, &n3, &n1, &n2, &n3, in, out, &sign); 
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
