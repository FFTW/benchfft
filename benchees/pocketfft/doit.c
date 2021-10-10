/* this program is in the public domain */

#include "bench-user.h"
#include "pocketfft.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "pocketfft")
BENCH_DOC("author", "Martin Reinecke")
BENCH_DOC("year", "2019")
BENCH_DOC("version", "1.0")
BENCH_DOC("language", "C")
BENCH_DOC("url", "https://gitlab.mpcdf.mpg.de/mtr/pocketfft")
BENCH_DOC("url-was-valid-on", "Fri Jul 23 23:06:24 ACST 2020")
BENCH_DOC("copyright", "3 clause BSDL")
END_BENCH_DOC


int can_do(struct problem *p)
{
     return (DOUBLE_PRECISION && p->rank == 1 && problem_in_place(p));
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_1d_fftpack(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_fftpack(p, in, -1.0);
}

static cfft_plan cplan;
static rfft_plan rplan;

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));

     if (p->kind == PROBLEM_COMPLEX) {
	 cplan = make_cfft_plan(p->n[0]);
     } else {
	 rplan = make_rfft_plan(p->n[0]);
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     int n = p->n[0];
     void *in = p->in;

     if (p->kind == PROBLEM_COMPLEX) {
	  if (p->sign == -1) {
	       for (i = 0; i < iter; ++i) {
		   cfft_forward(cplan, in, 1.0);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) {
		   cfft_backward(cplan, in, 1.0);
	       }
	  }
     } else {
	  if (p->sign == -1) {
	       for (i = 0; i < iter; ++i) {
		   rfft_forward(rplan, in, 1.0);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) {
		   rfft_backward(rplan, in, 1.0);
	       }
	  }
     }
}

void done(struct problem *p)
{
     if (p->kind == PROBLEM_COMPLEX)
	 destroy_cfft_plan(cplan);
     else
	 destroy_rfft_plan(rplan);
}
