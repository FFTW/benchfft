/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "cxml")
BENCH_DOC("package", "Compaq/Digital Extended Math Library (CXML/DXML)")
BENCH_DOC("notes", "3d transform")
BENCH_DOC("notes", "We benchmark using the complex array format.")
END_BENCH_DOC

#if defined(HAVE_DXMLDEF_H)
#  include <dxmldef.h>
#elif defined(HAVE_CXMLDEF_H)
#  include <cxmldef.h>
#endif

#define CFFT_INIT_3D F77_FUNC_(cfft_init_3d,CFFT_INIT_3D)
#define ZFFT_INIT_3D F77_FUNC_(zfft_init_3d,ZFFT_INIT_3D)
#define SFFT_INIT_3D F77_FUNC_(sfft_init_3d,SFFT_INIT_3D)
#define DFFT_INIT_3D F77_FUNC_(dfft_init_3d,DFFT_INIT_3D)
#define CFFT_APPLY_3D F77_FUNC_(cfft_apply_3d,CFFT_APPLY_3D)
#define ZFFT_APPLY_3D F77_FUNC_(zfft_apply_3d,ZFFT_APPLY_3D)
#define SFFT_APPLY_3D F77_FUNC_(sfft_apply_3d,SFFT_APPLY_3D)
#define DFFT_APPLY_3D F77_FUNC_(dfft_apply_3d,DFFT_APPLY_3D)
#define CFFT_EXIT_3D F77_FUNC_(cfft_exit_3d,CFFT_EXIT_3D)
#define ZFFT_EXIT_3D F77_FUNC_(zfft_exit_3d,ZFFT_EXIT_3D)
#define SFFT_EXIT_3D F77_FUNC_(sfft_exit_3d,SFFT_EXIT_3D)
#define DFFT_EXIT_3D F77_FUNC_(dfft_exit_3d,DFFT_EXIT_3D)

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
	       CFFT_INIT_3D(&ni, &nj, &nk, &fsc, &stride1);
	  else	       
	       ZFFT_INIT_3D(&ni, &nj, &nk, &fsz, &stride1);
     } else {
	  if (SINGLE_PRECISION) 
	       SFFT_INIT_3D(&ni, &nj, &nk, &fss, &stride1);
	  else	       
	       DFFT_INIT_3D(&ni, &nj, &nk, &fsd, &stride1);
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
		    CFFT_APPLY_3D("C", "C", dir, in, out, &lda, &ldb, &fsc, &stride,
				   &stride, &stride);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) {
		    ZFFT_APPLY_3D("C", "C", dir, in, out, &lda, &ldb, &fsz, &stride, 
				   &stride, &stride);
	       }
	  }
     } else {
	  int lda = p->n[2] + 2;
	  int ldb = p->n[1];
	  if (p->sign == -1) {
	       if (SINGLE_PRECISION) {
		    for (i = 0; i < iter; ++i) {
			 SFFT_APPLY_3D("R", "C", "F", in, out, &lda, &ldb, &fss, 
					&stride, &stride, &stride);
		    }
	       } else {
		    for (i = 0; i < iter; ++i) {
			 DFFT_APPLY_3D("R", "C", "F", in, out, &lda, &ldb, &fsd, 
					&stride, &stride, &stride);
		    }
	       }
	  } else {
	       if (SINGLE_PRECISION) {
		    for (i = 0; i < iter; ++i) {
			 SFFT_APPLY_3D("C", "R", "B", in, out, &lda, &ldb, &fss,
					&stride, &stride, &stride);
		    }
	       } else {
		    for (i = 0; i < iter; ++i) {
			 DFFT_APPLY_3D("C", "R", "B", in, out, &lda, &ldb, &fsd,
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
	       CFFT_EXIT_3D(&fsc); 
	  else          
	       ZFFT_EXIT_3D(&fsz); 
     } else { 
	  if (SINGLE_PRECISION)  
	       SFFT_EXIT_3D(&fss); 
	  else             
	       DFFT_EXIT_3D(&fsd); 
     } 
} 
 
