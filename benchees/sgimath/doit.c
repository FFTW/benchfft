/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "sgimath")
END_BENCH_DOC


int can_do(struct problem *p)
{
     return (p->rank == 1 && problem_in_place(p));
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_1d_fftpack(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_fftpack(p, in, -1.0);
}

#include <fft.h>

static void *WSAVE;

void setup(struct problem *p)
{
     int n;
 
     BENCH_ASSERT(can_do(p));
     n = p->n[0];
 
     if (p->kind == PROBLEM_COMPLEX) {
          WSAVE = bench_malloc((n + 15) * sizeof(bench_complex));
	  if (SINGLE_PRECISION) 
	       cfft1di(n, WSAVE);
	  else	       
	       zfft1di(n, WSAVE);
     } else {
          WSAVE = bench_malloc((n + 15) * sizeof(bench_real));
	  if (SINGLE_PRECISION) 
	       scfft1dui(n, WSAVE);
	  else	       
	       dzfft1dui(n, WSAVE);
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
	  if (SINGLE_PRECISION) {
	       for (i = 0; i < iter; ++i) {
		    cfft1d(sign, n, in, 1, wsave);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) {
		    zfft1d(sign, n, in, 1, wsave);
	       }
	  }
     } else {
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(WSAVE);
}
