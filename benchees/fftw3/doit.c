#include "bench-user.h"
#include <math.h>
#include <stdio.h>
#include <fftw3.h>
#include <string.h>

#define CONCAT(prefix, name) prefix ## name
#if defined(BENCHFFT_SINGLE)
#define FFTW(x) CONCAT(fftwf_, x)
#elif defined(BENCHFFT_LDOUBLE)
#define FFTW(x) CONCAT(fftwl_, x)
#else
#define FFTW(x) CONCAT(fftw_, x)
#endif

static const char *mkversion(void) { return FFTW(version); }
static const char *mkcc(void) { return FFTW(cc); }
static const char *mkcodelet_optim(void) { return FFTW(codelet_optim); }

BEGIN_BENCH_DOC
BENCH_DOC("name", "fftw3")
BENCH_DOCF("version", mkversion)
BENCH_DOCF("fftw-cc", mkcc)
BENCH_DOCF("fftw-codelet-optim", mkcodelet_optim)
END_BENCH_DOC 

FFTW(plan) the_plan = 0;
unsigned the_flags = 0;

void useropt(const char *arg)
{
     if (!strcmp(arg, "patient")) the_flags |= FFTW_PATIENT;
     else if (!strcmp(arg, "estimate")) the_flags |= FFTW_ESTIMATE;
     else if (!strcmp(arg, "exhaustive")) the_flags |= FFTW_EXHAUSTIVE;
     else if (!strcmp(arg, "unaligned")) the_flags |= FFTW_UNALIGNED;

     else fprintf(stderr, "unknown user option: %s.  Ignoring.\n", arg);
}

int can_do(struct problem *p)
{
     UNUSED(p);
     return 1;
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

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
 
     if (p->kind == PROBLEM_COMPLEX) {
	  the_plan = FFTW(plan_dft)(
	       p->rank, p->n,
	       p->in, p->out,
	       p->sign, the_flags);
     } else {
	  if (p->sign == -1) {
	       the_plan = FFTW(plan_dft_r2c)(
		    p->rank, p->n,
		    p->in, p->out,
		    the_flags);
	  } else {
	       the_plan = FFTW(plan_dft_c2r)(
		    p->rank, p->n,
		    p->in, p->out,
		    the_flags);
	  }
     }
     BENCH_ASSERT(the_plan);
}

void doit(int iter, struct problem *p)
{
     int i;
     FFTW(plan) plan = the_plan;

     UNUSED(p);

     for (i = 0; i < iter; ++i) 
	  FFTW(execute)(plan);
}

void done(struct problem *p)
{
     UNUSED(p);
     FFTW(destroy_plan)(the_plan);
}
