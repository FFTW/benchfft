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
     copy_h2c_1d_unpacked(p, out, -1.0);
}


void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_unpacked(p, in, -1.0);
}


void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, 1);
}

#ifdef HAVE_MKL_FFTC_LN_H
#include <mkl_fftc_ln.h>
#endif
#ifdef HAVE_MKL_FFT_H
#include <mkl_fft.h>
#endif

#define ZFFT1D zfft1d_
#define CFFT1D cfft1d_
#define DZFFT1D dzfft1d_
#define ZDFFT1D zdfft1d_
#define SCFFT1D scfft1d_
#define CSFFT1D csfft1d_

extern void ZFFT1D(), CFFT1D(), DZFFT1D(), ZDFFT1D(), SCFFT1D(), CSFFT1D();
static void *WSAVE;

void setup(struct problem *p)
{
     int n, zero = 0;
 
     BENCH_ASSERT(can_do(p));
     n = p->n[0];
 
     if (p->kind == PROBLEM_COMPLEX) {
	  /*
	   * example code says that wsave consists of 3 * n
	   * locations, but the code dumps core for n == 4
	   */
          WSAVE = bench_malloc((3 * n + 4) * sizeof(bench_real));
	  if (SINGLE_PRECISION) 
	       CFFT1D(p->in, &n, &zero, WSAVE);
	  else	       
	       ZFFT1D(p->in, &n, &zero, WSAVE);
     } else {
          WSAVE = bench_malloc((4 * n) * sizeof(bench_real));

	  if (p->sign == -1) {
	       if (SINGLE_PRECISION) 
		    SCFFT1D(p->in, &n, &zero, WSAVE);
	       else	       
		    DZFFT1D(p->in, &n, &zero, WSAVE);
	  } else {
	       if (SINGLE_PRECISION) 
		    CSFFT1D(p->in, &n, &zero, WSAVE);
	       else	       
		    ZDFFT1D(p->in, &n, &zero, WSAVE);
	  }
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     int n = p->n[0];
     void *wsave = WSAVE;
     int sign = p->sign;

     if (p->kind == PROBLEM_COMPLEX) {
	  bench_real *pin = p->in;

	  if (SINGLE_PRECISION) {
	       for (i = 0; i < iter; ++i) {
		    CFFT1D(pin, &n, &sign, wsave);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) {
		    ZFFT1D(pin, &n, &sign, wsave);
	       }
	  }
     } else {
	  bench_real *pin = p->in;

	  if (sign == -1) {
	       if (SINGLE_PRECISION) 
		    for (i = 0; i < iter; ++i) {
			 SCFFT1D(pin, &n, &sign, WSAVE);
		    }
	       else	       
		    for (i = 0; i < iter; ++i) {
			 DZFFT1D(pin, &n, &sign, WSAVE);
		    }
	  } else {
	       if (SINGLE_PRECISION) 
		    for (i = 0; i < iter; ++i) {
			 CSFFT1D(pin, &n, &sign, WSAVE);
		    }
	       else	       
		    for (i = 0; i < iter; ++i) {
			 ZDFFT1D(pin, &n, &sign, WSAVE);
		    }
	  }

     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(WSAVE);
}
