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
     copy_h2c_1d_unpacked_ri(p, out, -1.0);
}


void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_unpacked_ri(p, in, -1.0);
}


void copy_c2c_from(struct problem *p, bench_complex *in)
{
     copy_c2ri(in, p->in, ((bench_real *)p->in) + p->n[0], p->n[0]);
}

void copy_c2c_to(struct problem *p, bench_complex *out)
{
     copy_ri2c(p->out, ((bench_real *)p->out) + p->n[0], out, p->n[0]);
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, 1);
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

	  /*
	   * example code says that wsave consists of 3 * n
	   * locations, but the code dumps core for n == 4
	   */
          WSAVE = bench_malloc((3 * n + 4) * sizeof(bench_real));
	  if (SINGLE_PRECISION) 
	       CFFT1DC(rin, iin, n, 0, WSAVE);
	  else	       
	       ZFFT1DC(rin, iin, n, 0, WSAVE);
     } else {
          WSAVE = bench_malloc((4 * n) * sizeof(bench_real));

	  if (p->sign == -1) {
	       if (SINGLE_PRECISION) 
		    SCFFT1DC(p->in, n, 0, WSAVE);
	       else	       
		    DZFFT1DC(p->in, n, 0, WSAVE);
	  } else {
	       if (SINGLE_PRECISION) 
		    CSFFT1DC(p->in, n, 0, WSAVE);
	       else	       
		    ZDFFT1DC(p->in, n, 0, WSAVE);
	  }
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     int n = p->n[0];
     void *wsave = WSAVE;

     if (p->kind == PROBLEM_COMPLEX) {
	  bench_real *rin = p->in;
	  bench_real *iin = rin + n;
	  int sign = p->sign;

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
	  bench_real *pin = p->in;
	  if (p->sign == -1) {
	       if (SINGLE_PRECISION) 
		    for (i = 0; i < iter; ++i) {
			 SCFFT1DC(pin, n, -1, WSAVE);
		    }
	       else	       
		    for (i = 0; i < iter; ++i) {
			 DZFFT1DC(pin, n, -1, WSAVE);
		    }
	  } else {
	       if (SINGLE_PRECISION) 
		    for (i = 0; i < iter; ++i) {
			 CSFFT1DC(pin, n, 1, WSAVE);
		    }
	       else	       
		    for (i = 0; i < iter; ++i) {
			 ZDFFT1DC(pin, n, 1, WSAVE);
		    }
	  }

     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(WSAVE);
}
