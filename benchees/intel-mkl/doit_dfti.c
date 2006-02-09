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

#define ERROR_CHECK							\
     if (!DftiErrorClass(status, DFTI_NO_ERROR)) {			\
	  if (verbose > 3)						\
	       printf("DFTI error: %s\n", DftiErrorMessage(status));	\
	  return 0;							\
     }

#define MAXRNK 7 /* or so the paper claims */
static int mkdescriptor(struct problem *p)
{
     long status;
     enum DFTI_CONFIG_VALUE domain;
     long pn[MAXRNK];
     int i;

     the_descriptor = 0;

     if (p->rank > MAXRNK)
	  return 0;

     for (i = 0; i < p->rank; ++i)
	  pn[i] = p->n[i]; /* convert int -> long */

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
					pn[0]);
     else 
	  status = DftiCreateDescriptor(&the_descriptor,
					PRECISION,
					domain,
					p->rank,
					pn);
     ERROR_CHECK;

     status = DftiSetValue(the_descriptor, DFTI_PLACEMENT, 
			   problem_in_place(p) ? 
			   DFTI_INPLACE : DFTI_NOT_INPLACE);
     ERROR_CHECK;

     status = DftiSetValue(the_descriptor,
			   DFTI_CONJUGATE_EVEN_STORAGE, 
			   DFTI_COMPLEX_COMPLEX);
     ERROR_CHECK;

     /* we are already benchmarking CCS in doit.c */
     status = DftiSetValue(the_descriptor, DFTI_PACKED_FORMAT, 
			   DFTI_CCE_FORMAT);
     ERROR_CHECK;


     if (p->kind == PROBLEM_REAL) {
	  /* these guys must be kidding */
	  long strides[MAXRNK+1];

	  strides[p->rank] = 1;
	  if (p->rank > 0) {
	       strides[p->rank - 1] = pn[p->rank - 1] / 2 + 1;

	       for (i = p->rank - 2; i > 0; --i) 
		    strides[i] = strides[i+1] * pn[i];

	       /* if this is not zero the thing doesn't work */
	       strides[0] = 0;
	  }
	  
	  status = DftiSetValue(the_descriptor, 
				(p->sign > 0 ? 
				 DFTI_INPUT_STRIDES : DFTI_OUTPUT_STRIDES), 
				strides);
	  ERROR_CHECK;
     }

#if 0
     /* this disappeared in mkl 8 */
     DftiSetValue(the_descriptor, DFTI_INITIALIZATION_EFFORT, DFTI_HIGH);
     ERROR_CHECK;
#endif

     status = DftiCommitDescriptor(the_descriptor);
     ERROR_CHECK;

     return 1;
}

int can_do(struct problem *p)
{
     if (PRECISION < 0) return 0;

     /* ask mkl whether it can do it or not */
     if (mkdescriptor(p)) {
	  if (the_descriptor)
	       DftiFreeDescriptor(&the_descriptor);
	  return 1;
     } else 
	  return 0;
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
     int status = mkdescriptor(p);

     BENCH_ASSERT(status);
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
