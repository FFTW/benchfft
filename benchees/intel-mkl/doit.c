/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

#if HAVE_MKL_BLAS_H
#include <mkl_blas.h>
#endif

static const char *mkvers(void)
{
     static char buf[160];
     MKLGetVersionString(buf, 160);
     return buf;
}


BEGIN_BENCH_DOC
BENCH_DOC("name", "intel-mkl")
BENCH_DOC("package", "Intel Math Kernel Library (MKL)")
BENCH_DOCF("version", mkvers)
BENCH_DOC("notes", "Complex data are stored in separate real/imag arrays.")
BENCH_DOC("notes", "backward transform is scaled")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (p->rank >= 1 && p->rank <= 2 && problem_power_of_two(p, 1));
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

#ifdef BENCHFFT_SINGLE
#define _CFFT1DC cfft1dc
#define _SCFFT1DC scfft1dc
#define _CSFFT1DC csfft1dc
#define _CFFT2DC cfft2dc
#define _SCFFT2DC scfft2dc
#define _CSFFT2DC csfft2dc
#else
#define _CFFT1DC zfft1dc
#define _SCFFT1DC dzfft1dc
#define _CSFFT1DC zdfft1dc
#define _CFFT2DC zfft2dc
#define _SCFFT2DC dzfft2dc
#define _CSFFT2DC zdfft2dc
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
		    _CFFT1DC(rin, iin, n, 0, WSAVE);
	       } else {
		    WSAVE = bench_malloc((4 * n) * sizeof(bench_real));

		    if (p->sign == -1) {
			 _SCFFT1DC(p->in, n, 0, WSAVE);
		    } else {
			 _CSFFT1DC(p->in, n, 0, WSAVE);
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

		    for (i = 0; i < iter; ++i) {
			 _CFFT1DC(rin, iin, n, sign, wsave);
		    }
	       } else {
		    bench_real *pin = p->in;
		    if (p->sign == -1) {
			 for (i = 0; i < iter; ++i) {
			      _SCFFT1DC(pin, n, -1, WSAVE);
			 }
		    } else {
			 for (i = 0; i < iter; ++i) {
			      _CSFFT1DC(pin, n, 1, WSAVE);
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
		    
		    for (i = 0; i < iter; ++i) {
			 _CFFT2DC(rpin, ipin, m, n, sign);
		    }
	       } else {
		    bench_real *pin = p->in;
		    if (sign == -1) {
			 for (i = 0; i < iter; ++i) {
			      _SCFFT2DC(pin, m, n);
			 } 
		    } else {
			 for (i = 0; i < iter; ++i) {
			      _CSFFT2DC(pin, m, n);
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
