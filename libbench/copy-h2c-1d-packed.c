#include "bench.h"

/* 
 * copy a 1D packed hermitian array into a 1D complex array
 * Works only for even N.
 */

void copy_h2c_1d_packed(struct problem *p, bench_complex *out, 
			bench_real sign_of_r2h_transform)
{
     unsigned int k, n;
     bench_real *pout = p->out;

     BENCH_ASSERT(p->rank == 1);
     BENCH_ASSERT(p->kind == PROBLEM_REAL);

     n = p->n[0];
     BENCH_ASSERT((n & 1) == 0); /* even n */

     if (n > 0) {
	  c_re(out[0]) = pout[0];
	  c_im(out[0]) = 0.0;
     }

     for (k = 1; k < n / 2; ++k) {
	  c_re(out[k]) = pout[2 * k];
	  c_im(out[k]) = - sign_of_r2h_transform * pout[2 * k + 1];
	  c_re(out[n - k]) = pout[2 * k];
	  c_im(out[n - k]) = sign_of_r2h_transform * pout[2 * k + 1];
     }

     if (n > 0) {
	  c_re(out[k]) = pout[1];
	  c_im(out[k]) = 0.0;
     }
}
