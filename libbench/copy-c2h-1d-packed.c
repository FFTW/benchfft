#include "bench.h"

/* 
 * copy a 1D complex array into a 1D hermitian array stored in
 * ``packed'' format:
 *
 *    a[2*k] = R[k], 0<=k<n/2
 *    a[2*k+1] = I[k], 0<k<n/2
 *    a[1] = R[n/2]
 * 
 * Works only for even N.
 */

void copy_c2h_1d_packed(struct problem *p, bench_complex *in,
			bench_real sign_of_r2h_transform)
{
     unsigned int k, n;
     bench_real *pin = p->in;

     BENCH_ASSERT(p->rank == 1);
     BENCH_ASSERT(p->kind == PROBLEM_REAL);

     n = p->n[0];
     BENCH_ASSERT((n & 1) == 0); /* even n */

     for (k = 0; k < n / 2; ++k) {
	  pin[2 * k] = c_re(in[k]);
	  pin[2 * k + 1] = -sign_of_r2h_transform * c_im(in[k]);
     }
     if (k > 0)
	  pin[1] = c_re(in[k]);
}
