/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "sgimath")
BENCH_DOC("package", "Silicon Graphics Scientific Mathematical Library (complib.sgimath)")
END_BENCH_DOC


int can_do(struct problem *p)
{
     return (p->rank == 1 && problem_in_place(p));
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_unpacked(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_unpacked(p, in, -1.0);
}

#include <fft.h>

static void *WSAVE;

#ifdef BENCHFFT_SINGLE
#define CFFT1DI cfft1di
#define CFFT1D cfft1d
#define SCFFT1DUI scfft1dui
#define SCFFT1DU scfft1du
#define CSFFT1DU csfft1du
#else
#define CFFT1DI zfft1di
#define CFFT1D zfft1d
#define SCFFT1DUI dzfft1dui
#define SCFFT1DU dzfft1du
#define CSFFT1DU zdfft1du
#endif


void setup(struct problem *p)
{
     int n;
 
     BENCH_ASSERT(can_do(p));
     n = p->n[0];
 
     if (p->kind == PROBLEM_COMPLEX) {
          WSAVE = bench_malloc((n + 15) * sizeof(bench_complex));
	  CFFT1DI(n, WSAVE);
     } else {
          WSAVE = bench_malloc((n + 15) * sizeof(bench_real));
	  SCFFT1DUI(n, WSAVE); /* works for both directions */
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     int n = p->n[0];
     void *in = p->in;
     int sign = p->sign;
     void *wsave = WSAVE;

     if (p->kind == PROBLEM_COMPLEX) {
	  for (i = 0; i < iter; ++i) {
	       CFFT1D(sign, n, in, 1, wsave);
	  }
     } else {
	  if (p->sign == -1) {
	       for (i = 0; i < iter; ++i) {
		    SCFFT1DU(sign, n, in, 1,wsave);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) {
		    CSFFT1DU(sign, n, in, 1,wsave);
	       }
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(WSAVE);
}
