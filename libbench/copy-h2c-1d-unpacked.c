#include "bench.h"

/* 
 * copy a 1D unpacked hermitian array into a 1D complex array
 * (A unpacked hermitian array consists of the first 1 + n/2
 * elements of a complex array)
 */

void copy_h2c_1d_unpacked(struct problem *p, bench_complex *out, 
			  bench_real sign_of_r2h_transform)
{
     unsigned int k, n;
     bench_complex *pout = p->out;

     BENCH_ASSERT(p->rank == 1);
     BENCH_ASSERT(p->kind == PROBLEM_REAL);

     n = p->n[0];

     if (n > 0) {
	  c_re(out[0]) = c_re(pout[0]);
	  c_im(out[0]) = - sign_of_r2h_transform * c_im(pout[0]);
     }

     for (k = 1; k < 1 + (n / 2); ++k) {
	  c_re(out[k]) = c_re(pout[k]);
	  c_im(out[k]) = - sign_of_r2h_transform * c_im(pout[k]);
	  c_re(out[n - k]) = c_re(pout[k]);
	  c_im(out[n - k]) = sign_of_r2h_transform * c_im(pout[k]);
     }
}
