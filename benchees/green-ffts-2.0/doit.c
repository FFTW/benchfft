/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "green")
BENCH_DOC("package", "FFTs for RISC 2.0")
BENCH_DOC("year", "1998")
BENCH_DOC("author", "John Green")
BENCH_DOC("email", "green_jt@vsdec.npt.nuwc.navy.mil")
BENCH_DOC("copyright",
	  "This code is public domain, do anything you want to with it.")
BENCH_DOC("language", "C")
BENCH_DOC("notes", "The backward transform is scaled")
BENCH_DOC("notes", "Code is in single precision, but double precision version is created as documented by the author: ``for double precision just use a global search and replace to change float to double in all source files.''")
BENCH_DOC("url", "http://hyperarchive.lcs.mit.edu/HyperArchive/Archive/dev/src/ffts-for-risc-2-c.hqx")
BENCH_DOC("url-was-valid-on", "Sun Jul 15 21:28:59 EDT 2001")
END_BENCH_DOC

#include "fftext.h"
#include "fft2d.h"

int can_do(struct problem *p)
{
     return (problem_power_of_two(p, 1) &&
	     ((p->kind == PROBLEM_COMPLEX && p->rank >= 1 && p->rank < 4)
	      ||
	      (p->kind == PROBLEM_REAL && p->rank >= 1 && p->rank < 3)));
}


void copy_h2c(struct problem *p, bench_complex *out)
{
     switch (p->rank) {
	 case 1:
	      copy_h2c_1d_packed(p, out, -1.0);
	      break;
	 case 2:
	 {
	      bench_real *pout = p->out;
	      unsigned int i, j;
	      unsigned int n0 = p->n[0], n1 = p->n[1];

	      /* THIS FORMAT IS INSANE */

	      /* first half of first column */
	      for (i = 0; i < n0 / 2; ++i) {
		   c_re(out[i * n1]) = pout[i * n1];
		   c_im(out[i * n1]) = pout[i * n1 + 1];
	      }
	      c_im(out[0]) = 0;

	      /* nyquist value */
	      c_re(out[i * n1]) = pout[1];
	      c_im(out[i * n1]) = 0;
	      ++i;

	      /* expand first column by hermitian symmetry */
	      for (; i < n0; ++i) {
		   c_re(out[i * n1]) = c_re(out[(n0 - i) * n1]);
		   c_im(out[i * n1]) = -c_im(out[(n0 - i) * n1]);
	      }

	      /* first half of nyquist column */
	      for (i = 0; i < n0 / 2; ++i) {
		   c_re(out[i * n1 + n1 / 2]) = pout[(i + n0 / 2) * n1];
		   c_im(out[i * n1 + n1 / 2]) = pout[(i + n0 / 2) * n1 + 1];
	      }
	      c_im(out[n1 / 2]) = 0;

	      /* nyquist value */
	      c_re(out[i * n1 + n1 / 2]) = pout[n0 / 2 * n1 + 1];
	      c_im(out[i * n1 + n1 / 2]) = 0;
	      ++i;

	      /* expand nyquist column by hermitian symmetry */
	      for (; i < n0; ++i) {
		   c_re(out[i * n1 + n1 / 2]) =
			c_re(out[(n0 - i) * n1 + n1 / 2]);
		   c_im(out[i * n1 + n1 / 2]) = 
			-c_im(out[(n0 - i) * n1 + n1 / 2]);
	      }

	      /* remaining columns in first half */
	      for (j = 1; j < n1 / 2; ++j) {
		   for (i = 0; i < n0; ++i) {
			c_re(out[i * n1 + j]) = pout[i * n1 + 2 * j];
			c_im(out[i * n1 + j]) = pout[i * n1 + 2 * j + 1];
		   }
	      }

	      /* expand remaining columns by symmetry */
	      for (j = n1 / 2 + 1; j < n1; ++j) {
		   c_re(out[j]) = c_re(out[n1 - j]);
		   c_im(out[j]) = -c_im(out[n1 - j]);

		   for (i = 1; i < n0; ++i) {
			c_re(out[i * n1 + j]) = 
			     c_re(out[(n0 - i) * n1 + (n1 - j)]);
			c_im(out[i * n1 + j]) =
			     -c_im(out[(n0 - i) * n1 + (n1 - j)]);
		   }
	      }
	      break;
	 }
	 default:
	      BENCH_ASSERT(/* can't happen */ 0);
	      break;
     }
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     switch (p->rank) {
	 case 1:
	      copy_c2h_1d_packed(p, in, -1.0);
	      break;
	 case 2:
	 {
	      bench_real *pin = p->in;
	      unsigned int i, j;
	      unsigned int n0 = p->n[0], n1 = p->n[1];

	      /* THIS FORMAT IS INSANE */

	      /* first half of first column */
	      for (i = 0; i < n0 / 2; ++i) {
		   pin[i * n1] = c_re(in[i * n1]);
		   pin[i * n1 + 1] = c_im(in[i * n1]);
	      }

	      /* nyquist value */
	      pin[1] = c_re(in[i * n1]);

	      /* first half of nyquist column */
	      for (i = 0; i < n0 / 2; ++i) {
		   pin[(i + n0 / 2) * n1] = c_re(in[i * n1 + n1 / 2]);
		   pin[(i + n0 / 2) * n1 + 1] = c_im(in[i * n1 + n1 / 2]);
	      }

	      /* nyquist value */
	      pin[n0 / 2 * n1 + 1] = c_re(in[i * n1 + n1 / 2]);

	      /* remaining columns in first half */
	      for (j = 1; j < n1 / 2; ++j) {
		   for (i = 0; i < n0; ++i) {
			pin[i * n1 + 2 * j] = c_re(in[i * n1 + j]);
			pin[i * n1 + 2 * j + 1] = c_im(in[i * n1 + j]);
		   }
	      }
	      break;
	 }
	 default:
	      BENCH_ASSERT(/* can't happen */ 0);
	      break;
     }
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, 1);
}

static int m0, m1, m2;
void setup(struct problem *p)
{
     switch (p->rank) {
	 case 1: 
	      m0 = log_2(p->n[0]);
	      fftInit(m0);
	      break;
	 case 2: 
	      m0 = log_2(p->n[0]);
	      m1 = log_2(p->n[1]);
	      fft2dInit(m0, m1);
	      break;
	 case 3:
	      m0 = log_2(p->n[0]);
	      m1 = log_2(p->n[1]);
	      m2 = log_2(p->n[2]);
	      fft3dInit(m0, m1, m2); 
	      break;
	 default:
	      BENCH_ASSERT(/* can't happen */ 0);
	      break;
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     void *in = p->in;

     if (p->kind == PROBLEM_COMPLEX) {
	  switch (p->rank) {
	      case 1: 
	      {
		   if (p->sign == -1) {
			for (i = 0; i < iter; ++i) {
			     ffts(in, m0, 1);
			}
		   } else {
			for (i = 0; i < iter; ++i) {
			     iffts(in, m0, 1);
			}
		   }
		   break;
	      }
	      case 2: 
	      {
		   if (p->sign == -1) {
			for (i = 0; i < iter; ++i) {
			     fft2d(in, m0, m1);
			}
		   } else {
			for (i = 0; i < iter; ++i) {
			     ifft2d(in, m0, m1);
			}
		   }
		   break;
	      }
	      case 3: 
	      {
		   if (p->sign == -1) {
			for (i = 0; i < iter; ++i) {
			     fft3d(in, m0, m1, m2);
			}
		   } else {
			for (i = 0; i < iter; ++i) {
			     ifft3d(in, m0, m1, m2);
			}
		   }
		   break;
	      }
	  }
     } else {
	  switch (p->rank) {
	      case 1: 
	      {
		   if (p->sign == -1) {
			for (i = 0; i < iter; ++i) {
			     rffts(in, m0, 1);
			}
		   } else {
			for (i = 0; i < iter; ++i) {
			     riffts(in, m0, 1);
			}
		   }
		   break;
	      }
	      case 2: 
	      {
		   if (p->sign == -1) {
			for (i = 0; i < iter; ++i) {
			     rfft2d(in, m0, m1);
			}
		   } else {
			for (i = 0; i < iter; ++i) {
			     rifft2d(in, m0, m1);
			}
		   }
		   break;
	      }
	  }
     }
}

void done(struct problem *p)
{
     switch (p->rank) {
	 case 1: 
	      fftFree();
	      break;
	 case 2: 
	      fft2dFree();
	      break;
	 case 3:
	      fft3dFree();
	      break;
	 default:
	      BENCH_ASSERT(/* can't happen */ 0);
	      break;
     }
}
