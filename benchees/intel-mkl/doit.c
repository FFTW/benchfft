/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", NAME)
BENCH_DOC("notes", "backward transform is scaled")
END_BENCH_DOC

extern int have_sse(void);
static unsigned int cpuid_edx(unsigned int op)
{
     unsigned int eax, ebx, ecx, edx;

     __asm__("cpuid"
	     : "=a" (eax), "=b" (ebx), "=c" (ecx), "=d" (edx)
	     : "a" (op));
     return edx;
}

int can_do(struct problem *p)
{
     return (RIGHT_PROCESSOR &&
	     p->rank >= 1 && p->rank <= 2 && problem_power_of_two(p, 1));
}

/*

The following routine comes from Intel's library examples.  I am glad
I don't have to use this screwed format.  The routine is (C) Intel and
it is distributed accorging to the following excerpt from the license:

<p class=MsoPlainText><b>Sample Source:</b> may include example interface or
application source code.<span style="mso-spacerun: yes">  </span>You may copy,
modify and compile the Sample Source and distribute it in your own products in
binary and source code form.<o:p></o:p></p> 
*/

static void dz_copy(bench_real *a, bench_complex *o, int m, int n)
{
     int i, j;
     int n2, m2, nn;

     n2 = n + 2;
     m2 = m + 2;

     c_re(o[0]) = a[0];
     c_im(o[0]) = 0.0;

     if (m > 1) {
	  c_re(o[(m / 2) * n]) = a[(m / 2) * n2];
	  c_im(o[(m / 2) * n]) = 0.0;
     }

     for (i = 1; i < m / 2; i++) {
	  c_re(o[i * n]) = a[i * n2];
	  c_im(o[i * n]) = a[(m / 2 + i + 1) * n2];
	  c_re(o[(m - i) * n]) = a[i * n2];
	  c_im(o[(m - i) * n]) = -a[(m / 2 + i + 1) * n2];
     }

     if (n > 1) {
	  nn = n / 2;
	  c_re(o[nn]) = a[nn];
	  c_im(o[nn]) = 0.0;

	  if (m > 1) {
	       c_re(o[(m / 2) * n + nn]) = a[(m / 2) * n2 + nn];
	       c_im(o[(m / 2) * n + nn]) = 0.0;
	  }

	  for (i = 1; i < m / 2; i++) {
	       c_re(o[i * n + nn]) = a[i * n2 + nn];
	       c_im(o[i * n + nn]) = a[(m / 2 + i + 1) * n2 + nn];
	       c_re(o[(m - i) * n + nn]) = a[i * n2 + nn];
	       c_im(o[(m - i) * n + nn]) = -a[(m / 2 + i + 1) * n2 + nn];
	  }
     }


     for (j = 1; j < n / 2; j++) {
	  c_re(o[j]) = a[j];
	  c_im(o[j]) = a[n / 2 + j + 1];
	  c_re(o[n - j]) = a[j];
	  c_im(o[n - j]) = -a[n / 2 + j + 1];
	  for (i = 1; i < m; i++) {
	       c_re(o[i * n + j]) = a[i * n2 + j];
	       c_im(o[i * n + j]) = a[i * n2 + n / 2 + j + 1];
	       c_re(o[(m - i) * n + n - j]) = a[i * n2 + j];
	       c_im(o[(m - i) * n + n - j]) = -a[i * n2 + n / 2 + j + 1];
	  }
     }

}

/* complex->hermitian 2d copy routine */
static void zd_copy(bench_real *a, bench_complex *o, int m, int n)
{
     int i, j;
     int n2, m2, nn;

     n2 = n + 2;
     m2 = m + 2;

     a[0] = c_re(o[0]);

     if (m > 1) {
	  a[(m / 2) * n2] = c_re(o[(m / 2) * n]);
     }

     for (i = 1; i < m / 2; i++) {
	  a[i * n2] = c_re(o[i * n]);
	  a[(m / 2 + i + 1) * n2] = c_im(o[i * n]);
     }

     if (n > 1) {
	  nn = n / 2;
	  a[nn] = c_re(o[nn]);

	  if (m > 1) {
	       a[(m / 2) * n2 + nn] = c_re(o[(m / 2) * n + nn]);
	  }

	  for (i = 1; i < m / 2; i++) {
	       a[i * n2 + nn] = c_re(o[i * n + nn]);
	       a[(m / 2 + i + 1) * n2 + nn] = c_im(o[i * n + nn]);
	  }
     }


     for (j = 1; j < n / 2; j++) {
	  a[j] = c_re(o[j]);
	  a[n / 2 + j + 1] = c_im(o[j]);
	  for (i = 1; i < m; i++) {
	       a[i * n2 + j] = c_re(o[i * n + j]);
	       a[i * n2 + n / 2 + j + 1] = c_im(o[i * n + j]);
	  }
     }
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     if (p->rank == 1) {
	  copy_h2c_1d_unpacked_ri(p, out, -1.0);
     } else {
	  dz_copy(p->out, out, p->n[0], p->n[1]);
     }
}


void copy_c2h(struct problem *p, bench_complex *in)
{
     if (p->rank == 1) {
	  copy_c2h_1d_unpacked_ri(p, in, -1.0);
     } else {
	  zd_copy(p->in, in, p->n[0], p->n[1]);
     }
}


void copy_r2c(struct problem *p, bench_complex *out)
{
     copy_r2c_unpacked(p, out);
}

void copy_c2r(struct problem *p, bench_complex *in)
{
     copy_c2r_unpacked(p, in);
}

void copy_c2c_from(struct problem *p, bench_complex *in)
{
     copy_c2ri(in, p->in, ((bench_real *) p->in) + p->size, p->size);
}

void copy_c2c_to(struct problem *p, bench_complex *out)
{
     copy_ri2c(p->out, ((bench_real *) p->out) + p->size, out, p->size);
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

     switch (p->rank) {
     case 1:
	  {

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
	       break;
	  }

     case 2:
	  /* nothing to do */
	  break;
     default:
	  BENCH_ASSERT(0);
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     int sign = p->sign;

     switch (p->rank) {
     case 1:{
	       void *wsave = WSAVE;
	       int n = p->n[0];
	       if (p->kind == PROBLEM_COMPLEX) {
		    bench_real *rin = p->in;
		    bench_real *iin = rin + n;

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
	       break;
	  }
     case 2:
	  {
	       int m = p->n[0];
	       int n = p->n[1];

	       if (p->kind == PROBLEM_COMPLEX) {
		    bench_real *rpin = p->in;
		    bench_real *ipin = rpin + p->size;

		    if (SINGLE_PRECISION) {
			 for (i = 0; i < iter; ++i) {
			      CFFT2DC(rpin, ipin, m, n, sign);
			 }
		    } else {
			 for (i = 0; i < iter; ++i) {
			      ZFFT2DC(rpin, ipin, m, n, sign);
			 }
		    }
	       } else {
		    bench_real *pin = p->in;
		    if (sign == -1) {
			 if (SINGLE_PRECISION)
			      for (i = 0; i < iter; ++i) {
				   SCFFT2DC(pin, m, n);
			      } 
			 else
			      for (i = 0; i < iter; ++i) {
				   DZFFT2DC(pin, m, n);
			      }
		    } else {
			 if (SINGLE_PRECISION)
			      for (i = 0; i < iter; ++i) {
				   CSFFT2DC(pin, m, n);
			      } 
			 else
			      for (i = 0; i < iter; ++i) {
				   ZDFFT2DC(pin, m, n);
			      }
		    }
	       }
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(WSAVE);
}
