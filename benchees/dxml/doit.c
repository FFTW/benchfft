/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "dxml")
BENCH_DOC("notes", "using complex array format in real transforms")
END_BENCH_DOC


#include <dxmldef.h>

int can_do(struct problem *p)
{
     return (p->rank == 1 && 
	     /* N must be even for real transforms */
	     (p->kind == PROBLEM_COMPLEX || !(p->n[0] & 1)));
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_1d_unpacked(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_unpacked(p, in, -1.0);
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, 1);
}

static DXML_C_FFT_STRUCTURE fsc;
static DXML_Z_FFT_STRUCTURE fsz;
static DXML_S_FFT_STRUCTURE fss;
static DXML_D_FFT_STRUCTURE fsd;

void setup(struct problem *p)
{
     int n;
     int stride1 = 1;

     BENCH_ASSERT(can_do(p));
     n = p->n[0];
 
     if (p->kind == PROBLEM_COMPLEX) {
	  if (SINGLE_PRECISION) 
	       cfft_init_(&n, &fsc, &stride1);
	  else	       
	       zfft_init_(&n, &fsz, &stride1);
     } else {
	  if (SINGLE_PRECISION) 
	       sfft_init_(&n, &fss, &stride1);
	  else	       
	       dfft_init_(&n, &fsd, &stride1);
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
	  if (SINGLE_PRECISION) {
	       for (i = 0; i < iter; ++i) {
		    cfft_apply_("C", "C", dir, in, out, &fsc, &stride);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) {
		    zfft_apply_("C", "C", dir, in, out, &fsz, &stride);
	       }
	  }
     } else {
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
     }
}

void done(struct problem *p) 
{ 
     UNUSED(p); 
     if (p->kind == PROBLEM_COMPLEX) { 
	  if (SINGLE_PRECISION)  
	       cfft_exit_(&fsc); 
	  else          
	       zfft_exit_(&fsz); 
     } else { 
	  if (SINGLE_PRECISION)  
	       sfft_exit_(&fss); 
	  else             
	       dfft_exit_(&fsd); 
     } 
} 
 
