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
#include "doit-doc.h"
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (sizeof(fftw_real) == sizeof(bench_real) && (p->rank >= 1));
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_unpacked(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_unpacked(p, in, -1.0);
}

void copy_r2c(struct problem *p, bench_complex *out)
{
     if (problem_in_place(p))
	  copy_r2c_unpacked(p, out);	  
     else
	  copy_r2c_packed(p, out);
}

void copy_c2r(struct problem *p, bench_complex *in)
{
     if (problem_in_place(p))
	  copy_c2r_unpacked(p, in);
     else
	  copy_c2r_packed(p, in);
}

static void *plan;

void setup(struct problem *p)
{
     int n;
 
     BENCH_ASSERT(can_do(p));
 
     if (p->kind == PROBLEM_COMPLEX) {
	  plan = fftwnd_create_plan_specific(
	       p->rank, p->n, p->sign,
	       (PLAN_FLAGS | (problem_in_place(p) ? FFTW_IN_PLACE : 0)),
	       p->in, 1,
	       p->out, 1);
     } else {
	  plan = rfftwnd_create_plan_specific(
	       p->rank, p->n, p->sign,
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
	       fftwnd_one(PLAN, in, out);
	  }
     } else {
	  if (p->sign == -1)
	       for (i = 0; i < iter; ++i) {
		    rfftwnd_one_real_to_complex(PLAN, in, out);
	       }
	  else	
	       for (i = 0; i < iter; ++i) {
		    rfftwnd_one_complex_to_real(PLAN, in, out);
	       }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     if (p->kind == PROBLEM_COMPLEX) {
	  fftwnd_destroy_plan(plan);
     } else {
	  rfftwnd_destroy_plan(plan);
     }
}
