/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "goedecker")
END_BENCH_DOC

#define SG_FFT_ F77_FUNC_(sg_fft, SG_FFT)
extern void SG_FFT_();
 
int can_do(struct problem *p)
{
     return (p->rank == 3 &&
	     DOUBLE_PRECISION &&
	     p->kind == PROBLEM_COMPLEX &&
	     p->n[0] <= 2048 && p->n[1] <= 2048 && p->n[2] <= 2048 &&
	     check_prime_factors(p->size, 5) &&
	     !problem_in_place(p));
}

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
     /* nothing to do */
}

void doit(int iter, struct problem *p)
{
     int i;
     int n1 = p->n[2];
     int n2 = p->n[1];
     int n3 = p->n[0];
     int fftcache = 16;
     void *in = p->in;
     void *out = p->out;
     double sign = p->sign;

     for (i = 0; i < iter; ++i) {
	  SG_FFT_(&fftcache, &n1, &n2, &n3, &n1, &n2, &n3, in, out, &sign); 
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
