#include "bench-user.h"
#include <math.h>
#include <stdio.h>
#include <mkl_dfti.h>

#ifdef BENCHFFT_SINGLE
#define PRECISION DFTI_SINGLE
#elif defined(BENCHFFT_LDOUBLE)
#define PRECISION (-1) /* can't do it */
#else
#define PRECISION DFTI_DOUBLE
#endif

#if HAVE_MKL_BLAS_H
#include <mkl_blas.h>
#endif

static const char *mkvers(void)
{
#if HAVE_MKLGETVERSIONSTRING
     static char buf[160];
     MKLGetVersionString(buf, 160);
     return buf;
#else
     return "unknown"
#endif
}

BEGIN_BENCH_DOC
BENCH_DOC("name", "mkl-dfti")
BENCH_DOC("package", "Intel Math Kernel Library (MKL), DFTI interface")
BENCH_DOCF("version", mkvers)
END_BENCH_DOC

DFTI_DESCRIPTOR *the_descriptor;

int can_do(struct problem *p)
{
     /* real transforms appear to be unimplemented in mkl6-beta */
     return (p->kind == PROBLEM_COMPLEX &&
	     sizeof(bench_real) <= sizeof(double));
}

void setup(struct problem *p)
{
     long status;
     enum DFTI_CONFIG_VALUE domain;

     BENCH_ASSERT(can_do(p));
 
     if (p->kind == PROBLEM_COMPLEX) {
	  domain = DFTI_COMPLEX;
     } else {
	  if (p->sign == -1)
	       domain = DFTI_REAL;
	  else
	       domain = DFTI_CONJUGATE_EVEN;
     }
     
     /* api nonsense */
     if (p->rank == 1) 
	  status = DftiCreateDescriptor(&the_descriptor,
					PRECISION,
					domain,
					1,
					p->n[0]);
     else 
	  status = DftiCreateDescriptor(&the_descriptor,
					PRECISION,
					domain,
					p->rank,
					p->n);

     if (status) {
	  printf("DFTI error: %s\n", DftiErrorMessage(status));
	  exit(1);
     }
     DftiSetValue(the_descriptor, DFTI_PLACEMENT, 
		  problem_in_place(p) ? DFTI_INPLACE : DFTI_NOT_INPLACE);
     DftiSetValue(the_descriptor, DFTI_INITIALIZATION_EFFORT, DFTI_HIGH);
		  
     DftiCommitDescriptor(the_descriptor);
}

void doit(int iter, struct problem *p)
{
     int i;
     DFTI_DESCRIPTOR *d = the_descriptor;

     if (p->in_place) {
	  if (p->sign == -1) {
	       for (i = 0; i < iter; ++i) 
		    DftiComputeForward(d, p->in);
	  } else {
	       for (i = 0; i < iter; ++i) 
		    DftiComputeBackward(d, p->in);
	  }
     } else {
	  if (p->sign == -1) {
	       for (i = 0; i < iter; ++i) 
		    DftiComputeForward(d, p->in, p->out);
	  } else {
	       for (i = 0; i < iter; ++i) 
		    DftiComputeBackward(d, p->in, p->out);
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     DftiFreeDescriptor(&the_descriptor);
}
