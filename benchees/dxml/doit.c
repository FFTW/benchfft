/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "cxml")
BENCH_DOC("package", "Compaq/Digital Extended Math Library (CXML/DXML)")
BENCH_DOC("notes", "1d transform")
BENCH_DOC("notes", "We benchmark using the complex array format.")
END_BENCH_DOC

#if defined(HAVE_DXMLDEF_H)
#  include <dxmldef.h>
#elif defined(HAVE_CXMLDEF_H)
#  include <cxmldef.h>
#endif

int can_do(struct problem *p)
{
     return (p->rank == 1 && 
	     /* N must be even for real transforms */
	     (p->kind == PROBLEM_COMPLEX || !(p->n[0] & 1)));
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_unpacked(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_unpacked(p, in, -1.0);
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, 1);
}

#define CFFT_INIT F77_FUNC_(cfft_init,CFFT_INIT)
#define ZFFT_INIT F77_FUNC_(zfft_init,ZFFT_INIT)
#define SFFT_INIT F77_FUNC_(sfft_init,SFFT_INIT)
#define DFFT_INIT F77_FUNC_(dfft_init,DFFT_INIT)
#define CFFT_APPLY F77_FUNC_(cfft_apply,CFFT_APPLY)
#define ZFFT_APPLY F77_FUNC_(zfft_apply,ZFFT_APPLY)
#define SFFT_APPLY F77_FUNC_(sfft_apply,SFFT_APPLY)
#define DFFT_APPLY F77_FUNC_(dfft_apply,DFFT_APPLY)
#define CFFT_EXIT F77_FUNC_(cfft_exit,CFFT_EXIT)
#define ZFFT_EXIT F77_FUNC_(zfft_exit,ZFFT_EXIT)
#define SFFT_EXIT F77_FUNC_(sfft_exit,SFFT_EXIT)
#define DFFT_EXIT F77_FUNC_(dfft_exit,DFFT_EXIT)

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
	       CFFT_INIT(&n, &fsc, &stride1);
	  else	       
	       ZFFT_INIT(&n, &fsz, &stride1);
     } else {
	  if (SINGLE_PRECISION) 
	       SFFT_INIT(&n, &fss, &stride1);
	  else	       
	       DFFT_INIT(&n, &fsd, &stride1);
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
	  if (SINGLE_PRECISION) {
	       for (i = 0; i < iter; ++i) {
		    CFFT_APPLY("C", "C", dir, in, out, &fsc, &stride);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) {
		    ZFFT_APPLY("C", "C", dir, in, out, &fsz, &stride);
	       }
	  }
     } else {
	  if (p->sign == -1) {
	       if (SINGLE_PRECISION) {
		    for (i = 0; i < iter; ++i) {
			 SFFT_APPLY("R", "C", "F", in, out, &fss, &stride);
		    }
	       } else {
		    for (i = 0; i < iter; ++i) {
			 DFFT_APPLY("R", "C", "F", in, out, &fsd, &stride);
		    }
	       }
	  } else {
	       if (SINGLE_PRECISION) {
		    for (i = 0; i < iter; ++i) {
			 SFFT_APPLY("C", "R", "B", in, out, &fss, &stride);
		    }
	       } else {
		    for (i = 0; i < iter; ++i) {
			 DFFT_APPLY("C", "R", "B", in, out, &fsd, &stride);
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
	       CFFT_EXIT(&fsc); 
	  else          
	       ZFFT_EXIT(&fsz); 
     } else { 
	  if (SINGLE_PRECISION)  
	       SFFT_EXIT(&fss); 
	  else             
	       DFFT_EXIT(&fsd); 
     } 
} 
 
