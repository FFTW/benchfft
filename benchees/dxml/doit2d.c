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
	     (p->kind == PROBLEM_COMPLEX || !(p->n[1] & 1)));
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, 1);
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
     copy_r2c_unpacked(p, out);	  
}

void copy_c2r(struct problem *p, bench_complex *in)
{
     copy_c2r_unpacked(p, in);
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
     void *in = p->in;
     void *out = p->out;
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
	  int lda = p->n[1] + 2; /* n[1] is even */
	  if (p->sign == -1) {
	       if (SINGLE_PRECISION) {
		    for (i = 0; i < iter; ++i) {
			 sfft_apply_2d_("R", "C", "F", in, out, &lda, &fss, 
					&stride, &stride);
		    }
	       } else {
		    for (i = 0; i < iter; ++i) {
			 dfft_apply_2d_("R", "C", "F", in, out, &lda, &fsd, 
					&stride, &stride);
		    }
	       }
	  } else {
	       if (SINGLE_PRECISION) {
		    for (i = 0; i < iter; ++i) {
			 sfft_apply_2d_("C", "R", "B", in, out, &lda, &fss, 
					&stride, &stride);
		    }
	       } else {
		    for (i = 0; i < iter; ++i) {
			 dfft_apply_2d_("C", "R", "B", in, out, &lda, &fsd,
					&stride, &stride);
		    }
	       }
	  }
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
 