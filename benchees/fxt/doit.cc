/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

#include "fxt.h"

#include "doit-doc.h"

#ifndef RANK_OK
#  define RANK_OK(r) ((r) == 1)
#endif

extern "C"
int can_do(struct problem *p)
{
     return (DOUBLE_PRECISION &&
	     p->kind == PROBLEM_COMPLEX &&
	     RANK_OK(p->rank) &&
	     problem_power_of_two(p, 1));
}

Complex *x = (Complex *) 0;
int m, m1, N, n0, n1;
ulong *md;

extern "C"
void setup(struct problem *p)
{
     int i, n = p->n[0];
#ifndef NO_Complex
     x = new Complex[n];
     for (i = 0; i < n; ++i) {
	  x[i] = Complex(0.0, 0.0);
     }
#endif
     md = new ulong[p->rank];
     N = 1;
     for (i = 0; i < p->rank; ++i) {
       md[i] = log_2(p->n[p->rank - 1 - i]);
       N *= p->n[i];
     }
     if (p->rank > 0) {
       n0 = p->n[0];
       m = md[p->rank - 1 - 0];
     }
     if (p->rank > 1) {
       n1 = p->n[1];
       m1 = md[p->rank - 1 - 1];
     }
}

#ifndef NO_Complex
extern "C"
void problem_ccopy_from(struct problem *p, bench_complex *in)
{
     int i;
     Complex *out = x;

     for (i = 0; i < N; ++i) {
	  out[i] = Complex(c_re(in[i]), c_im(in[i]));
     }
}

extern "C"
void problem_ccopy_to(struct problem *p, bench_complex *out)
{
     int i;
     Complex *in = x;

     for (i = 0; i < N; ++i) {
	  c_re(out[i]) = real(in[i]);
	  c_im(out[i]) = imag(in[i]);
     }
}
#else /* NO_Complex */
void problem_ccopy_from(struct problem *p, bench_complex *in)
{
     bench_real *x = (bench_real *) p->in;
     copy_c2ri(in, x, x + N, N);
}

void problem_ccopy_to(struct problem *p, bench_complex *out)
{
     bench_real *x = (bench_real *) p->in;
     copy_ri2c(x, x + N, out, N);
}
#endif /* NO_Complex */

extern "C"
void doit(int iter, struct problem *p)
{
     int i;
     int is = p->sign;
#ifdef NO_Complex
     int rank= p->rank;
     double *inr = (double *) p->in;
     double *ini = inr + N;
#endif

     for (i = 0; i < iter; ++i) {
	  DOIT_FFT;
     }
}

extern "C"
void done(struct problem *p)
{
     UNUSED(p);
     delete md;
#ifndef NO_Complex
     delete x;
#endif
}
