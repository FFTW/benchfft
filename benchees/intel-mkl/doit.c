/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", NAME)
BENCH_DOC("notes", "backward transform is scaled")
END_BENCH_DOC


int can_do(struct problem *p)
{
     return (p->rank == 1 && problem_power_of_two(p, 1));
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_1d_fftpack(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_fftpack(p, in, -1.0);
}

void problem_ccopy_from(struct problem *p, bench_complex *in)
{
     BENCH_ASSERT(p->kind == PROBLEM_COMPLEX && p->rank == 1);
     copy_c2ri(in, p->in, ((bench_real *)p->in) + p->n[0],
	       p->n[0]);
}

void problem_ccopy_to(struct problem *p, bench_complex *out)
{
     BENCH_ASSERT(p->kind == PROBLEM_COMPLEX && p->rank == 1);
     copy_ri2c(p->out, ((bench_real *)p->out) + p->n[0], out,
	       p->n[0]);

     if (p->sign == 1) { 	  /* undo the scaling */
	  bench_complex x;
	  c_re(x) = p->size;
	  c_im(x) = 0;

	  cascale(out, p->size, x);
     }
}

#ifdef HAVE_MKL_FFT_H
#include <mkl_fft.h>
#endif

#ifdef HAVE_MKL_FFTC_LN_H
#include <mkl_fftc_ln.h>
#endif

static void *WSAVE;

void setup(struct problem *p)
{
     int n;
 
     BENCH_ASSERT(can_do(p));
     n = p->n[0];
 
     if (p->kind == PROBLEM_COMPLEX) {
	  bench_real *rin = p->in;
	  bench_real *iin = rin + n;

          WSAVE = bench_malloc((3 * n) * sizeof(bench_real));
	  if (SINGLE_PRECISION) 
	       CFFT1DC(rin, iin, n, 0, WSAVE);
	  else	       
	       ZFFT1DC(rin, iin, n, 0, WSAVE);
     } else {
#if 0
          WSAVE = bench_malloc((n + 15) * sizeof(bench_real));
	  if (SINGLE_PRECISION) 
	       scfft1dui(n, WSAVE);
	  else	       
	       dzfft1dui(n, WSAVE);
#endif
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     int n = p->n[0];
     bench_real *rin = p->in;
     bench_real *iin = rin + n;
     int sign = p->sign;
     void *wsave = WSAVE;

     if (p->kind == PROBLEM_COMPLEX) {
	  if (SINGLE_PRECISION) {
	       for (i = 0; i < iter; ++i) {
		    CFFT1DC(rin, iin, n, sign, wsave);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) {
		    ZFFT1DC(rin, iin, n, sign, wsave);
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
