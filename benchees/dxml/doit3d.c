/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "dxml")
BENCH_DOC("package", "Compaq/Digital Extended Math Library (CXML/DXML)")
BENCH_DOC("notes", "3d transform")
BENCH_DOC("notes", "We benchmark using the complex array format.")
END_BENCH_DOC


#include <dxmldef.h>

int can_do(struct problem *p)
{
     return (p->rank == 3 && 
	     /* N must be even for real transforms */
	     (p->kind == PROBLEM_COMPLEX || !(p->n[2] & 1)));
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

static DXML_C_FFT_STRUCTURE_3D fsc;
static DXML_Z_FFT_STRUCTURE_3D fsz;
static DXML_S_FFT_STRUCTURE_3D fss;
static DXML_D_FFT_STRUCTURE_3D fsd;

void setup(struct problem *p)
{
     int ni;
     int nj;
     int nk;
     int stride1 = 1;

     BENCH_ASSERT(can_do(p));
     ni = p->n[2];
     nj = p->n[1];
     nk = p->n[0];
 
     if (p->kind == PROBLEM_COMPLEX) {
	  if (SINGLE_PRECISION) 
	       cfft_init_3d_(&ni, &nj, &nk, &fsc, &stride1);
	  else	       
	       zfft_init_3d_(&ni, &nj, &nk, &fsz, &stride1);
     } else {
	  if (SINGLE_PRECISION) 
	       sfft_init_3d_(&ni, &nj, &nk, &fss, &stride1);
	  else	       
	       dfft_init_3d_(&ni, &nj, &nk, &fsd, &stride1);
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
	  int lda = p->n[2];
	  int ldb = p->n[1];
	  if (SINGLE_PRECISION) {
	       for (i = 0; i < iter; ++i) {
		    cfft_apply_3d_("C", "C", dir, in, out, &lda, &ldb, &fsc, &stride,
				   &stride, &stride);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) {
		    zfft_apply_3d_("C", "C", dir, in, out, &lda, &ldb, &fsz, &stride, 
				   &stride, &stride);
	       }
	  }
     } else {
	  int lda = p->n[2] + 2;
	  int ldb = p->n[1];
	  if (p->sign == -1) {
	       if (SINGLE_PRECISION) {
		    for (i = 0; i < iter; ++i) {
			 sfft_apply_3d_("R", "C", "F", in, out, &lda, &ldb, &fss, 
					&stride, &stride, &stride);
		    }
	       } else {
		    for (i = 0; i < iter; ++i) {
			 dfft_apply_3d_("R", "C", "F", in, out, &lda, &ldb, &fsd, 
					&stride, &stride, &stride);
		    }
	       }
	  } else {
	       if (SINGLE_PRECISION) {
		    for (i = 0; i < iter; ++i) {
			 sfft_apply_3d_("C", "R", "B", in, out, &lda, &ldb, &fss,
					&stride, &stride, &stride);
		    }
	       } else {
		    for (i = 0; i < iter; ++i) {
			 dfft_apply_3d_("C", "R", "B", in, out, &lda, &ldb, &fsd,
					&stride, &stride, &stride);
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
	       cfft_exit_3d_(&fsc); 
	  else          
	       zfft_exit_3d_(&fsz); 
     } else { 
	  if (SINGLE_PRECISION)  
	       sfft_exit_3d_(&fss); 
	  else             
	       dfft_exit_3d_(&fsd); 
     } 
} 
 
