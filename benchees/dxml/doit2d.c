/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "dxml")
BENCH_DOC("notes", "using complex array format")
END_BENCH_DOC


#include <dxmldef.h>

int can_do(struct problem *p)
{
     return (p->rank == 2 && 
	     /* N must be even for real transforms */
	     (p->kind == PROBLEM_COMPLEX || !(p->n[0] & 1)));
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, 1);
}

static DXML_C_FFT_STRUCTURE_2D fsc;
static DXML_Z_FFT_STRUCTURE_2D fsz;
static DXML_S_FFT_STRUCTURE_2D fss;
static DXML_D_FFT_STRUCTURE_2D fsd;

void setup(struct problem *p)
{
     int ni;
     int nj;
     int stride1 = 1;

     BENCH_ASSERT(can_do(p));
     ni = p->n[1];
     nj = p->n[0];
 
     if (p->kind == PROBLEM_COMPLEX) {
	  if (SINGLE_PRECISION) 
	       cfft_init_2d_(&ni, &nj, &fsc, &stride1);
	  else	       
	       zfft_init_2d_(&ni, &nj, &fsz, &stride1);
     } else {
	  if (SINGLE_PRECISION) 
	       sfft_init_2d_(&ni, &nj, &fss, &stride1);
	  else	       
	       dfft_init_2d_(&ni, &nj, &fsd, &stride1);
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     int n = p->n[0];
     void *in = p->in;
     void *out = p->out;
     int sign = p->sign;
     int stride = 1;

     if (p->kind == PROBLEM_COMPLEX) {
	  char *dir = p->sign == -1 ? "F" : "B";
	  int lda = p->n[1];
	  if (SINGLE_PRECISION) {
	       for (i = 0; i < iter; ++i) {
		    cfft_apply_2d_("C", "C", dir, in, out, &lda, &fsc, &stride,
				   &stride);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) {
		    zfft_apply_2d_("C", "C", dir, in, out, &lda, &fsz, &stride, 
				   &stride);
	       }
	  }
     } else {
#if 0
	  if (p->sign == -1) {
	       if (SINGLE_PRECISION) {
		    for (i = 0; i < iter; ++i) {
			 sfft_apply_("R", "C", "F", in, out, &fss, &stride);
		    }
	       } else {
		    for (i = 0; i < iter; ++i) {
			 dfft_apply_("R", "C", "F", in, out, &fsd, &stride);
		    }
	       }
	  } else {
	       if (SINGLE_PRECISION) {
		    for (i = 0; i < iter; ++i) {
			 sfft_apply_("C", "R", "B", in, out, &fss, &stride);
		    }
	       } else {
		    for (i = 0; i < iter; ++i) {
			 dfft_apply_("C", "R", "B", in, out, &fsd, &stride);
		    }
	       }
	  }
#endif
     }
}

void done(struct problem *p) 
{ 
     UNUSED(p); 
     if (p->kind == PROBLEM_COMPLEX) { 
	  if (SINGLE_PRECISION)  
	       cfft_exit_2d_(&fsc); 
	  else          
	       zfft_exit_2d_(&fsz); 
     } else { 
	  if (SINGLE_PRECISION)  
	       sfft_exit_2d_(&fss); 
	  else             
	       dfft_exit_2d_(&fsd); 
     } 
} 
 
