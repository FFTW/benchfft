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
BENCH_DOC("name", "intel-mkl-f")
BENCH_DOC("package", "Intel Math Kernel Library (MKL)")
BENCH_DOCF("version", mkvers)
BENCH_DOC("notes", "backward transform is scaled")
END_BENCH_DOC 

int can_do(struct problem *p)
{
     return (p->rank >= 1 && p->rank <= 2 && problem_power_of_two(p, 1));
}


/*

The following routine comes from Intel's library examples.  I am glad
I don't have to use this screwed format.  The routine was translated
by f2c and massaged by hand.  The routine is (C) Intel and it is
distributed accorging to the following excerpt from the license:


<p class=MsoPlainText><b>Sample Source:</b> may include example interface or
application source code.<span style="mso-spacerun: yes">  </span>You may copy,
modify and compile the Sample Source and distribute it in your own products in
binary and source code form.<o:p></o:p></p>

 */

static void dzmcopy(bench_real *ra, int rlda, bench_complex *ca,
		    int clda, int m, int n)
{
     /* System generated locals */
     int ra_dim1, ra_offset, ca_dim1, ca_offset, i1, i2, i3, i4, i5;
     bench_real d1;

     /* Local variables */
     int i, j, k, mm;

/* --- Input Arguments */
/* --- Output Arguments */
/* --- Local variables */
     /* Parameter adjustments */
     ra_dim1 = rlda;
     ra_offset = 1 + ra_dim1 * 1;
     ra -= ra_offset;
     ca_dim1 = clda;
     ca_offset = 1 + ca_dim1 * 1;
     ca -= ca_offset;

     /* Function Body */
     i1 = ca_dim1 + 1;
     i2 = ra_dim1 + 1;
     ca[i1].re = ra[i2], ca[i1].im = 0.;
     if (n > 1) {
	  i1 = (n / 2 + 1) * ca_dim1 + 1;
	  i2 = (n + 1) * ra_dim1 + 1;
	  ca[i1].re = ra[i2], ca[i1].im = 0.;
     }
     k = 3;
     i1 = n / 2;
     for (j = 2; j <= i1; ++j) {
	  i2 = j * ca_dim1 + 1;
	  i3 = k * ra_dim1 + 1;
	  i4 = (k + 1) * ra_dim1 + 1;
	  ca[i2].re = ra[i3], ca[i2].im = ra[i4];
	  i2 = (n + 2 - j) * ca_dim1 + 1;
	  i3 = k * ra_dim1 + 1;
	  d1 = -ra[(k + 1) * ra_dim1 + 1];
	  ca[i2].re = ra[i3], ca[i2].im = d1;
	  k += 2;
     }
     if (m > 1) {
	  mm = m / 2 + 1;
	  i1 = mm + ca_dim1;
	  i2 = m + 1 + ra_dim1;
	  ca[i1].re = ra[i2], ca[i1].im = 0;
	  if (n > 1) {
	       i1 = mm + (n / 2 + 1) * ca_dim1;
	       i2 = m + 1 + (n + 1) * ra_dim1;
	       ca[i1].re = ra[i2], ca[i1].im = 0;
	  }
	  k = 3;
	  i1 = n / 2;
	  for (j = 2; j <= i1; ++j) {
	       i2 = mm + j * ca_dim1;
	       i3 = m + 1 + k * ra_dim1;
	       i4 = m + 1 + (k + 1) * ra_dim1;
	       ca[i2].re = ra[i3], ca[i2].im = ra[i4];
	       i2 = mm + (n + 2 - j) * ca_dim1;
	       i3 = m + 1 + k * ra_dim1;
	       d1 = -ra[m + 1 + (k + 1) * ra_dim1];
	       ca[i2].re = ra[i3], ca[i2].im = d1;
	       k += 2;
	  }
     }
     k = 3;
     i1 = m / 2;
     for (i = 2; i <= i1; ++i) {
	  i2 = i + ca_dim1;
	  i3 = k + ra_dim1;
	  i4 = k + 1 + ra_dim1;
	  ca[i2].re = ra[i3], ca[i2].im = ra[i4];
	  i2 = m + 2 - i + ca_dim1;
	  i3 = k + ra_dim1;
	  d1 = -ra[k + 1 + ra_dim1];
	  ca[i2].re = ra[i3], ca[i2].im = d1;
	  i2 = n;
	  for (j = 2; j <= i2; ++j) {
	       i3 = i + j * ca_dim1;
	       i4 = k + j * ra_dim1;
	       i5 = k + 1 + j * ra_dim1;
	       ca[i3].re = ra[i4], ca[i3].im = ra[i5];
	       i3 = m + 2 - i + (n + 2 - j) * ca_dim1;
	       i4 = k + j * ra_dim1;
	       d1 = -ra[k + 1 + j * ra_dim1];
	       ca[i3].re = ra[i4], ca[i3].im = d1;
	  }
	  k += 2;
     }
}				/* dzmcopy_ */


static void zdmcopy(bench_real *ra, int rlda, bench_complex *ca,
		    int clda, int m, int n)
{
     /* System generated locals */
     int ra_dim1, ra_offset, ca_dim1, ca_offset, i1, i2, i3, i4, i5;

     /* Local variables */
     int i, j, k, mm;

     /* Parameter adjustments */
     ra_dim1 = rlda;
     ra_offset = 1 + ra_dim1 * 1;
     ra -= ra_offset;
     ca_dim1 = clda;
     ca_offset = 1 + ca_dim1 * 1;
     ca -= ca_offset;

     /* Function Body */
     i1 = ca_dim1 + 1;
     i2 = ra_dim1 + 1;
     ra[i2] = ca[i1].re;
     if (n > 1) {
	  i1 = (n / 2 + 1) * ca_dim1 + 1;
	  i2 = (n + 1) * ra_dim1 + 1;
	  ra[i2] = ca[i1].re;
     }
     k = 3;
     i1 = n / 2;
     for (j = 2; j <= i1; ++j) {
	  i2 = j * ca_dim1 + 1;
	  i3 = k * ra_dim1 + 1;
	  i4 = (k + 1) * ra_dim1 + 1;
	  ra[i3] = ca[i2].re, ra[i4] = ca[i2].im;
	  k += 2;
     }
     if (m > 1) {
	  mm = m / 2 + 1;
	  i1 = mm + ca_dim1;
	  i2 = m + 1 + ra_dim1;
	  ra[i2] = ca[i1].re;
	  if (n > 1) {
	       i1 = mm + (n / 2 + 1) * ca_dim1;
	       i2 = m + 1 + (n + 1) * ra_dim1;
	       ra[i2] = ca[i1].re;
	  }
	  k = 3;
	  i1 = n / 2;
	  for (j = 2; j <= i1; ++j) {
	       i2 = mm + j * ca_dim1;
	       i3 = m + 1 + k * ra_dim1;
	       i4 = m + 1 + (k + 1) * ra_dim1;
	       ra[i3] = ca[i2].re, ra[i4] = ca[i2].im;
	       k += 2;
	  }
     }
     k = 3;
     i1 = m / 2;
     for (i = 2; i <= i1; ++i) {
	  i2 = i + ca_dim1;
	  i3 = k + ra_dim1;
	  i4 = k + 1 + ra_dim1;
	  ra[i3] = ca[i2].re, ra[i4] = ca[i2].im;
	  i2 = n;
	  for (j = 2; j <= i2; ++j) {
	       i3 = i + j * ca_dim1;
	       i4 = k + j * ra_dim1;
	       i5 = k + 1 + j * ra_dim1;
	       ra[i4] = ca[i3].re, ra[i5] = ca[i3].im;
	  }
	  k += 2;
     }
}


void copy_h2c(struct problem *p, bench_complex *out)
{
     if (p->rank == 1) {
	  copy_h2c_unpacked(p, out, -1.0);
     } else {
	  dzmcopy(p->out, p->n[1] + 2, out, p->n[1], p->n[1], p->n[0]);
     }
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     if (p->rank == 1) {
	  copy_c2h_unpacked(p, in, -1.0);
     } else {
	  zdmcopy(p->in, p->n[1] + 2, in, p->n[1], p->n[1], p->n[0]);
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

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, 1);
}

#ifdef HAVE_MKL_FFT_H
#include <mkl_fft.h>
#endif

#define ZFFT1D zfft1d_
#define CFFT1D cfft1d_
#define DZFFT1D dzfft1d_
#define ZDFFT1D zdfft1d_
#define SCFFT1D scfft1d_
#define CSFFT1D csfft1d_

#define ZFFT2D zfft2d_
#define CFFT2D cfft2d_
#define DZFFT2D dzfft2d_
#define ZDFFT2D zdfft2d_
#define SCFFT2D scfft2d_
#define CSFFT2D csfft2d_

extern void ZFFT1D(), CFFT1D(), DZFFT1D(), ZDFFT1D(), SCFFT1D(), CSFFT1D();
extern void ZFFT2D(), CFFT2D(), DZFFT2D(), ZDFFT2D(), SCFFT2D(), CSFFT2D();
static void *WSAVE;

void setup(struct problem *p)
{
     int n, zero = 0;

     BENCH_ASSERT(can_do(p));

     switch (p->rank) {
     case 1:
	  {
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
	       int n = p->n[0];
	       void *wsave = WSAVE;

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
	       break;
	  }
     case 2:{
	       int n = p->n[0];
	       int m = p->n[1];
	       bench_real *pin = p->in;

	       if (p->kind == PROBLEM_COMPLEX) {
		    if (SINGLE_PRECISION) {
			 for (i = 0; i < iter; ++i) {
			      CFFT2D(pin, &m, &n, &sign);
			 }
		    } else {
			 for (i = 0; i < iter; ++i) {
			      ZFFT2D(pin, &m, &n, &sign);
			 }
		    }
	       } else {
		    if (sign == -1) {
			 if (SINGLE_PRECISION)
			      for (i = 0; i < iter; ++i) {
				   SCFFT2D(pin, &m, &n);
			      }
			 else
			      for (i = 0; i < iter; ++i) {
				   DZFFT2D(pin, &m, &n);
			      }
		    } else {
			 if (SINGLE_PRECISION)
			      for (i = 0; i < iter; ++i) {
				   CSFFT2D(pin, &m, &n);
			 } else
			      for (i = 0; i < iter; ++i) {
				   ZDFFT2D(pin, &m, &n);
			      }
		    }
	       }
	  }

     }
}

void done(struct problem *p)
{
     UNUSED(p);
     if (p->rank == 1)
	  bench_free(WSAVE);
}
