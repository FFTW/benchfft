/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "sgimath")
BENCH_DOC("package", "Silicon Graphics Scientific Mathematical Library (complib.sgimath)")
END_BENCH_DOC


int can_do(struct problem *p)
{
     return (p->rank == 2 && problem_in_place(p));
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

#include <fft.h>


#ifdef BENCHFFT_SINGLE
#define CFFT2DI cfft2di
#define CFFT2D cfft2d
#define SCFFT2DUI scfft2dui
#define SCFFT2DU scfft2du
#define CSFFT2DU csfft2du
#else
#define CFFT2DI zfft2di
#define CFFT2D zfft2d
#define SCFFT2DUI dzfft2dui
#define SCFFT2DU dzfft2du
#define CSFFT2DU zdfft2du
#endif

static void *WSAVE;

void setup(struct problem *p)
{
     int n1, n2;
 
     BENCH_ASSERT(can_do(p));
     n1 = p->n[1];
     n2 = p->n[0];
 
     if (p->kind == PROBLEM_COMPLEX) {
          WSAVE = bench_malloc((n1 + n2 + 30) * sizeof(bench_complex));
	  CFFT2DI(n1, n2, WSAVE);
     } else {
          WSAVE = bench_malloc((n1 + 2 * n2 + 45) * sizeof(bench_real));
	  SCFFT2DUI(n1, n2, WSAVE); /* works for both directions */
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     int n1 = p->n[1];
     int n2 = p->n[0];
     void *in = p->in;
     int sign = p->sign;
     void *wsave = WSAVE;

     if (p->kind == PROBLEM_COMPLEX) {
	  for (i = 0; i < iter; ++i) {
	       CFFT2D(sign, n1, n2, in, n1, wsave);
	  }
     } else {
          int lda = 2 * (1 + n1 / 2);
	  if (p->sign == -1) {
	       for (i = 0; i < iter; ++i) {
		    SCFFT2DU(sign, n1, n2, in, lda, wsave);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) {
		    CSFFT2DU(sign, n1, n2, in, lda, wsave);
	       }
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(WSAVE);
}
