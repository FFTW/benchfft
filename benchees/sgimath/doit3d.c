/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "sgimath")
BENCH_DOC("package", "Silicon Graphics Scientific Mathematical Library (complib.sgimath)")
END_BENCH_DOC


int can_do(struct problem *p)
{
     return (p->rank == 3 && problem_in_place(p));
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

static void *WSAVE;


#ifdef BENCHFFT_SINGLE
#define CFFT3DI cfft3di
#define CFFT3D cfft3d
#define SCFFT3DUI scfft3dui
#define SCFFT3DU scfft3du
#define CSFFT3DU csfft3du
#else
#define CFFT3DI zfft3di
#define CFFT3D zfft3d
#define SCFFT3DUI dzfft3dui
#define SCFFT3DU dzfft3du
#define CSFFT3DU zdfft3du
#endif

void setup(struct problem *p)
{
     int n1, n2, n3;
 
     BENCH_ASSERT(can_do(p));
     n1 = p->n[2];
     n2 = p->n[1];
     n3 = p->n[0];
 
     if (p->kind == PROBLEM_COMPLEX) {
          WSAVE = bench_malloc((n1 + n2 + n3 + 45) * sizeof(bench_complex));
	  CFFT3DI(n1, n2, n3, WSAVE);
     } else {
          WSAVE = bench_malloc((n1 + 2 * n2 + 2 * n3 + 75)
			       * sizeof(bench_real));
	  SCFFT3DUI(n1, n2, n3, WSAVE); /* works for both directions */
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     int n1 = p->n[2];
     int n2 = p->n[1];
     int n3 = p->n[0];
     void *in = p->in;
     int sign = p->sign;
     void *wsave = WSAVE;

     if (p->kind == PROBLEM_COMPLEX) {
	  for (i = 0; i < iter; ++i) {
	       CFFT3D(sign, n1, n2, n3, in, n1, n2, wsave);
	  }
     } else {
          int lda = 2 * (1 + n1 / 2);
	  if (p->sign == -1) {
	       for (i = 0; i < iter; ++i) {
		    SCFFT3DU(sign, n1, n2, n3, in, lda, n2, wsave);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) {
		    CSFFT3DU(sign, n1, n2, n3, in, lda, n2, wsave);
	       }
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(WSAVE);
}
