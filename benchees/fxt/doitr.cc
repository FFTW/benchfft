/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

#include "fxt.h"

#include "doit-doc.h"

extern "C"
int can_do(struct problem *p)
{
     return (DOUBLE_PRECISION &&
	     p->kind == PROBLEM_REAL &&
	     p->rank == 1 &&
	     p->n[0] > 1 &&
	     problem_power_of_two(p, 1));
}

Complex *x = (Complex *) 0;
int m;

extern "C"
void setup(struct problem *p)
{
     m = log_2(p->n[0]);
}


#ifdef PACKED /* crazy format, not matching documentation */
extern "C"
void copy_h2c(struct problem *p, bench_complex *out)
{
  double *pout = (double *) p->out;
  int i, n = p->n[0];

  BENCH_ASSERT(n % 2 == 0);
  c_re(out[0]) = pout[0];
  c_re(out[n/2]) = pout[1];
  c_im(out[0]) = c_im(out[n/2]) = 0.0;
  for (i = 1; i < n - i; ++i) {
    c_re(out[i]) = c_re(out[n-i]) = pout[2*i];
    c_im(out[n-i]) = -(c_im(out[i]) = pout[n-2*i+1]);
  }
  // documented format: copy_h2c_1d_packed(p, out, -1.0);
}

extern "C"
void copy_c2h(struct problem *p, bench_complex *in)
{
  double *pin = (double *) p->in;
  int i, n = p->n[0];

  /* Grrr, c2r *same* sign as r2c */
  BENCH_ASSERT(n % 2 == 0);
  pin[0] = c_re(in[0]);
  pin[1] = c_re(in[n/2]);
  for (i = 1; i < n - i; ++i) {
    pin[2*i] = c_re(in[i]);
    pin[n-2*i+1] = -c_im(in[i]);
  }
  //documented format: copy_c2h_1d_packed(p, in, +1.0);
}
#endif /* PACKED */

#ifdef HALFCOMPLEX
extern "C"
void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_1d_halfcomplex(p, out, -1.0);
}

extern "C"
void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_halfcomplex(p, in, +1.0); /* Grrr, c2r *same* sign as r2c */
}
#endif /* HALFCOMPLEX */

extern "C"
void doit(int iter, struct problem *p)
{
     int i;
     double *in = (double *) p->in;
     
     if (p->sign == -1)
       for (i = 0; i < iter; ++i) {
	 DOIT_FFT;
       }
     else
       for (i = 0; i < iter; ++i) {
	 DOIT_IFFT;
       }
}

extern "C"
void done(struct problem *p)
{
     UNUSED(p);
}
