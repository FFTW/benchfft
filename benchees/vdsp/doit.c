/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

#include <vDSP.h>

/* The vec_malloc and vec_free functions are recommended in the Apple
   vDSP docs, but I can't find them anywhere; let's hope that malloc
   produces sufficiently-aligned arrays? */
#if !defined(HAVE_VEC_MALLOC) || !defined(HAVE_VEC_FREE)
#  define vec_malloc bench_malloc
#  define vec_free bench_free
#endif

BEGIN_BENCH_DOC
BENCH_DOC("name", "vdsp")
BENCH_DOC("author", "Apple Computer, Inc.")
BENCH_DOC("copyright", "Bundled with MacOS 9.1 and later versions.")
BENCH_DOC("notes",
	  "Part of Apple's vecLib, this library is designed for use with "
	  "the AltiVec SIMD instructions on a PowerPC G4.")
BENCH_DOC("notes", "The forward real transform is scaled by 2.")
BENCH_DOC("notes", "Complex data are stored in separate real/imag arrays.")
BENCH_DOC("notes", "Real data are split into even/odd arrays.")
BENCH_DOC("url", "http://developer.apple.com/techpubs/macosx/CoreTechnologies/vDSP/vDSP.html")
BENCH_DOC("url-was-valid-on", "Tue Jul 17 23:01:08 EDT 2001")
BENCH_DOC("bibitem",
	  "R. Crandall and J. Klivington, Supercomputer-style FFT library "
	  "for Apple G4, Apple Technical Report (Jan. 2000).")
END_BENCH_DOC

#define CONCAT(prefix, name) prefix ## name

#ifdef BENCHFFT_SINGLE
#define MANGLE(name) name
typedef DSPSplitComplex splitcomplex;
#else
#define MANGLE(name) CONCAT(name, D)
typedef DSPDoubleSplitComplex splitcomplex;
#endif

static int n_ok(struct problem *p, unsigned int n)
{
     /* According to Apple's docs, the code only vectorizes for n
	between 2^2 and 2^20, inclusive, and other sizes allegedly use
	scalar code.  For real transforms, the limits start at 2^3.
	With MacOS X's vecLib, however, other sizes seem to crash or
	(worse) hang (sigh). */
     if (p->kind == PROBLEM_COMPLEX && (n < (1<<2) || n > (1<<20)))
	  return 0;
     if (p->kind == PROBLEM_REAL && (n < (1<<3) || n > (1<<20)))
	  return 0;
     return (power_of_two(n));
}

int can_do(struct problem *p)
{
     return (n_ok(p, p->n[0]) && (p->rank < 2 || n_ok(p, p->n[1])) &&
	     ((p->kind == PROBLEM_COMPLEX && p->rank >= 1 && p->rank < 3)  ||
	      (p->kind == PROBLEM_REAL && p->rank >= 1 && p->rank < 3 &&
	       p->n[0] > 1 && (p->rank < 2 || p->n[0] > 1))));
}

static int imax2(int a, int b) { return (a > b ? a : b); }
static int imin2(int a, int b) { return (a < b ? a : b); }

MANGLE(FFTSetup) fftsetup;
splitcomplex ins, outs, buf;
int m0, m1;
int n0, n1;

void setup(struct problem *p)
{
     n1 = p->rank == 2 ? p->n[1] : 1;
     fftsetup = MANGLE(create_fftsetup)(imax2(m0 = log_2(n0 = p->n[0]), 
					      m1 = log_2(n1)),
					2);

     /* Use Apple vec_malloc for 16-byte alignment */
     ins.realp = (bench_real*) vec_malloc(p->size * sizeof(bench_real));
     ins.imagp = (bench_real*) vec_malloc(p->size * sizeof(bench_real));
     if (!problem_in_place(p)) {
	  outs.realp = (bench_real*) vec_malloc(p->size * sizeof(bench_real));
	  outs.imagp = (bench_real*) vec_malloc(p->size * sizeof(bench_real));
     }
     else {
	  outs.realp = ins.realp;
	  outs.imagp = ins.imagp;
	  buf.realp = (bench_real*) vec_malloc(imin2(4 * p->size, 16384)
					       * sizeof(bench_real));
	  buf.imagp = (bench_real*) vec_malloc(imin2(4 * p->size, 16384)
					       * sizeof(bench_real));
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     FFTDirection dir = p->sign < 0 ? FFT_FORWARD : FFT_INVERSE;

     if (p->kind == PROBLEM_COMPLEX) {
	  if (problem_in_place(p)) {
	       switch (p->rank) {
		   case 1: 
			for (i = 0; i < iter; ++i) {
			     MANGLE(fft_zipt)(fftsetup, &ins, 1, &buf, m0, dir);
			}
			break;
		   case 2: 
			for (i = 0; i < iter; ++i) {
			     MANGLE(fft2d_zipt)(fftsetup, &ins, 1, 0, &buf, 
						m1, m0, dir);
			}
			break;
	       }
	  }
	  else {
	       switch (p->rank) {
		   case 1: 
			for (i = 0; i < iter; ++i) {
			     MANGLE(fft_zop)(fftsetup, &ins, 1, &outs, 1, m0, dir);
			}
			break;
		   case 2: 
			for (i = 0; i < iter; ++i) {
			     MANGLE(fft2d_zop)(fftsetup, &ins, 1, 0, &outs, 1, 0,
					       m1, m0, dir);
			}
			break;
	       }
	  }
     }
     else { /* PROBLEM_REAL */
	  if (problem_in_place(p)) {
	       switch (p->rank) {
		   case 1: 
			for (i = 0; i < iter; ++i) {
			     MANGLE(fft_zript)(fftsetup, &ins, 1, &buf, m0, dir);
			}
			break;
		   case 2: 
			for (i = 0; i < iter; ++i) {
			     MANGLE(fft2d_zript)(fftsetup, &ins, 1, 0, &buf, 
						 m1, m0, dir);
			}
			break;
	       }
	  }
	  else {
	       switch (p->rank) {
		   case 1: 
			for (i = 0; i < iter; ++i) {
			     MANGLE(fft_zrop)(fftsetup, &ins, 1, &outs, 1, m0, dir);
			}
			break;
		   case 2: 
			for (i = 0; i < iter; ++i) {
			     MANGLE(fft2d_zrop)(fftsetup, &ins, 1, 0, &outs, 1, 0,
						m1, m0, dir);
			}
			break;
	       }
	  }
     }
}

void done(struct problem *p)
{
     if (!problem_in_place(p)) {
	  vec_free(outs.imagp);
	  vec_free(outs.realp);
     }
     else {
	  vec_free(buf.imagp);
	  vec_free(buf.realp);
     }
     vec_free(ins.imagp);
     vec_free(ins.realp);
     MANGLE(destroy_fftsetup)(fftsetup);
}

void copy_c2c_from(struct problem *p, bench_complex *in)
{
     copy_c2ri(in, (bench_real*) ins.realp, (bench_real*) ins.imagp, p->size);
}

void copy_c2c_to(struct problem *p, bench_complex *out)
{
     copy_ri2c((bench_real*) outs.realp, (bench_real*) outs.imagp, 
	       out, p->size);
}

void copy_r2c(struct problem *p, bench_complex *out)
{
     unsigned int i;
     unsigned int n = p->size;

     for (i = 0; 2*i+1 < n; ++i) {
          c_re(out[2*i]) = outs.realp[i];
          c_im(out[2*i]) = 0.0;
          c_re(out[2*i+1]) = outs.imagp[i];
          c_im(out[2*i+1]) = 0.0;
     }
     if (2*i < n) {
          c_re(out[2*i]) = outs.realp[i];
          c_im(out[2*i]) = 0.0;
     }
}

void copy_c2r(struct problem *p, bench_complex *in)
{
     unsigned int i;
     unsigned int n = p->size;

     for (i = 0; 2*i+1 < n; ++i) {
          ins.realp[i] = c_re(in[2*i]);
	  ins.imagp[i] = c_re(in[2*i+1]);
     }
     if (2*i < n) {
          ins.realp[i] = c_re(in[2*i]);
     }
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     bench_complex *poutc = (bench_complex *) p->out;
     switch (p->rank) {
	 case 1:
	      copy_ri2c((bench_real*) outs.realp, (bench_real*) outs.imagp,
			poutc, p->size / 2);
	      copy_h2c_1d_packed(p, out, -1.0);
	      break;
	 case 2:
	 {
	      unsigned int i, j;
              unsigned int n0 = p->n[0], n1 = p->n[1];

	      for (i = 0; i < n0; ++i) {
		   for (j = 0; j < n1/2; ++j) {
			c_re(poutc[i*(n1/2+1)+j]) = outs.realp[i*(n1/2)+j];
			c_im(poutc[i*(n1/2+1)+j]) = outs.imagp[i*(n1/2)+j];
		   }
	      }

	      /* unpack Nyquist column */
	      for (i = 0; i < n0; ++i)
		   outs.realp[i] = c_im(poutc[i*(n1/2+1)]);
	      copy_h2c_1d_packed_strided(n0, (bench_real*) outs.realp, 1,
					 poutc + n1/2, n1/2+1, -1.0);

	      /* unpack DC column */
	      for (i = 0; i < n0; ++i)
		   outs.realp[i] = c_re(poutc[i*(n1/2+1)]);
	      copy_h2c_1d_packed_strided(n0, (bench_real*) outs.realp, 1,
					 poutc, n1/2+1, -1.0);
	      
	      copy_h2c_unpacked(p, out, -1.0);

	      break;
	 }
	 default:
	      BENCH_ASSERT(/* can't happen */ 0);
	      break;
     }
     {
	  bench_complex x;
	  c_re(x) = 0.5;
          c_im(x) = 0;
          cascale(out, p->size, x);
     }
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     bench_complex *pinc = (bench_complex *) p->in;
     switch (p->rank) {
	 case 1:
	      copy_c2h_1d_packed(p, in, -1.0);
	      copy_c2ri(pinc, (bench_real*) ins.realp, (bench_real*) ins.imagp,
			p->size / 2);
	      break;
	 case 2:
	 {
	      unsigned int i, j;
              unsigned int n0 = p->n[0], n1 = p->n[1];

	      copy_c2h_unpacked(p, in, -1.0);

	      /* unpack DC column (in place) */
	      copy_c2h_1d_packed_strided(n0, (bench_real*) ins.realp, 1,
					 pinc, n1/2+1, -1.0);
	      for (i = 0; i < n0; ++i)
		   c_re(pinc[i*(n1/2+1)]) = ins.realp[i];
	      
	      /* unpack Nyquist column */
	      copy_c2h_1d_packed_strided(n0, (bench_real*) ins.realp, 1,
					 pinc + n1/2, n1/2+1, -1.0);
	      for (i = 0; i < n0; ++i)
		   c_im(pinc[i*(n1/2+1)]) = ins.realp[i];

	      for (i = 0; i < n0; ++i) {
		   for (j = 0; j < n1/2; ++j) {
			ins.realp[i*(n1/2)+j] = c_re(pinc[i*(n1/2+1)+j]);
			ins.imagp[i*(n1/2)+j] = c_im(pinc[i*(n1/2+1)+j]);
		   }
	      }

	      break;
	 }
	 default:
	      BENCH_ASSERT(/* can't happen */ 0);
	      break;
     }
}


