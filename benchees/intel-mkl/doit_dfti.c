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
     static char buf[160];
     MKLGetVersionString(buf, 160);
     return buf;
}

BEGIN_BENCH_DOC
BENCH_DOC("name", "intel-mkl-dfti")
BENCH_DOC("package", "Intel Math Kernel Library (MKL), DFTI interface")
BENCH_DOCF("version", mkvers)
BENCH_DOC("notes", "Using the default CCS storage for hermitian data.")
END_BENCH_DOC

DFTI_DESCRIPTOR *the_descriptor;

static int mkdescriptor(struct problem *p)
{
     long status;
     enum DFTI_CONFIG_VALUE domain;

     if (p->kind == PROBLEM_COMPLEX) {
	  domain = DFTI_COMPLEX;
     } else {
	  domain = DFTI_REAL;
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

     if (status && verbose > 3) 
	  printf("DFTI error: %s\n", DftiErrorMessage(status));
     
     return (!status);
}

int can_do(struct problem *p)
{
     if (PRECISION < 0) return 0;

     /* ask mkl whether it can do it or not */
     if (mkdescriptor(p)) {
	  DftiFreeDescriptor(&the_descriptor);
	  return 1;
     } else 
	  return 0;
}


void copy_h2c(struct problem *p, bench_complex *out)
{
     if (p->rank == 1) {
	  copy_h2c_unpacked(p, out, -1.0);
     } else {
	  BENCH_ASSERT(0  /* not implemented in mkl 6.0 */);
     }
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     if (p->rank == 1) {
	  copy_c2h_unpacked(p, in, -1.0);
     } else {
	  BENCH_ASSERT(0  /* not implemented in mkl 6.0 */);
     }
}

void copy_r2c(struct problem *p, bench_complex *out)
{
     copy_r2c_unpacked(p, out);
}

void copy_c2r(struct problem *p, bench_complex *in)
{
     copy_c2r_unpacked(p, in);
}

void setup(struct problem *p)
{
     int status = mkdescriptor(p);

     BENCH_ASSERT(status);
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
