/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "djbfft-0.76")
BENCH_DOC("author", "D. J. Bernstein")
BENCH_DOC("year", "1999")
BENCH_DOC("email", "djb@pobox.com")
BENCH_DOC("url", "http://pobox.com/~djb/djbfft.html")
BENCH_DOC("url-was-valid-on", "Thu Jul 12 20:26:24 EDT 2001")
BENCH_DOC("language", "C")
BENCH_DOC("notes", "the forward transform has sign +1")
BENCH_DOC("notes", "the output of the forward transform is out of order")
BENCH_DOC("djbfft-compiled-by", DJBFFT_CC)
BENCH_DOC("conf-opt", CONF_OPT)
BENCH_DOC("auto-opt", AUTO_OPT)
END_BENCH_DOC

#include "fftc4.h"
#include "fftc8.h"
#include "fftr4.h"
#include "fftr8.h"
#include "fftfreq.h"

static void (*fft)();

int can_do(struct problem *p)
{
     return p->rank == 1 && problem_power_of_two(p, 1) && p->n[0] <= 8192;
}

void copy_c2c_from(struct problem *p, bench_complex *in)
{
     bench_complex *pin = p->in;

     if (p->sign == -1)
	  /* in-order input for forward transforms */
	  cacopy(in, p->in, p->n[0]);
     else {
	  unsigned int i;
	  for (i = 0; i < p->n[0]; ++i) {
	       pin[i] = in[fftfreq_c(i, p->n[0])];
	  }
     }

     /* conjugate the input because of opposite sign convention */
     {
	  unsigned int i;
	  for (i = 0; i < p->n[0]; ++i) {
	       c_im(pin[i]) = -c_im(pin[i]);
	  }
     }
}

void copy_c2c_to(struct problem *p, bench_complex *out)
{
     bench_complex *pout = p->out;

     /* conjugate the output because of opposite sign convention */
     {
	  unsigned int i;
	  for (i = 0; i < p->n[0]; ++i) {
	       c_im(pout[i]) = -c_im(pout[i]);
	  }
     }

     if (p->sign == -1) {
	  unsigned int i;
	  for (i = 0; i < p->n[0]; ++i) {
	       out[fftfreq_c(i, p->n[0])] = pout[i];
	  }
     } else {
	  /* in-order output for forward transforms */
	  cacopy(p->out, out, p->n[0]);
     }
}

/* interleaved input and out of order output.  What a f*cking idiot */
void copy_c2r(struct problem *p, bench_complex *in)
{
     unsigned int i;
     unsigned int n = p->size;
     bench_real *pin = p->in;

     for (i = 0; i < n / 2; ++i) {
	  pin[2 * i] = c_re(in[i]);
	  pin[2 * i + 1] = c_re(in[i + n / 2]);
     }
}

void copy_r2c(struct problem *p, bench_complex *out)
{
     unsigned int i;
     unsigned int n = p->size;
     bench_real *pout = p->out;

     for (i = 0; i < n / 2; ++i) {
	  c_re(out[i]) = pout[2 * i]; 
	  c_im(out[i]) = 0;
	  c_re(out[i + n / 2]) = pout[2 * i + 1];
	  c_im(out[i + n / 2]) = 0;
     }
}

void copy_h2c(struct problem *p, bench_complex *out) 
{
     unsigned int k, n;
     bench_real *pout = p->out;

     BENCH_ASSERT(p->rank == 1);
     BENCH_ASSERT(p->kind == PROBLEM_REAL);

     n = p->n[0];

     if (n > 0) {
	  c_re(out[0]) = pout[0];
	  c_im(out[0]) = 0;

	  c_re(out[n / 2]) = pout[1];
	  c_im(out[n / 2]) = 0;
     }

     for (k = 2; k < n; k += 2) {
          int f = fftfreq_r(k, n);
	  c_re(out[f]) = pout[k];
	  c_im(out[f]) = -pout[k + 1];
	  c_re(out[n - f]) = pout[k];
	  c_im(out[n - f]) = pout[k + 1];
     }
}

void copy_c2h(struct problem *p, bench_complex *in) 
{
     unsigned int k, n;
     bench_real *pin = p->in;

     BENCH_ASSERT(p->rank == 1);
     BENCH_ASSERT(p->kind == PROBLEM_REAL);

     n = p->n[0];

     if (n > 0) {
	  pin[0] = c_re(in[0]);
	  pin[1] = c_re(in[n / 2]);
     }

     for (k = 2; k < n; k += 2) {
          int f = fftfreq_r(k, n);
	  pin[k] = 2.0 * c_re(in[f]);
	  pin[k + 1] = - 2.0 * c_im(in[f]);
     }
}

#define CASE(p, var, fft)			\
switch (p->n[0]) {				\
    case 2: var =fft ## 2; break;		\
    case 4: var = fft ## 4; break;		\
    case 8: var = fft ## 8; break;		\
    case 16: var = fft ## 16; break;		\
    case 32: var = fft ## 32; break;		\
    case 64: var = fft ## 64; break;		\
    case 128: var = fft ## 128; break;		\
    case 256: var = fft ## 256; break;		\
    case 512: var = fft ## 512; break;		\
    case 1024: var = fft ## 1024; break;	\
    case 2048: var = fft ## 2048; break;	\
    case 4096: var = fft ## 4096; break;	\
    case 8192: var = fft ## 8192; break;	\
}

void setup(struct problem *p)
{
     if (p->kind == PROBLEM_COMPLEX) {
	  if (SINGLE_PRECISION) {
	       if (p->sign == -1) {
		    CASE(p, fft, fftc4_);
	       } else {
		    CASE(p, fft, fftc4_un);
	       }
	  } else {
	       if (p->sign == -1) {
		    CASE(p, fft, fftc8_);
	       } else {
		    CASE(p, fft, fftc8_un);
	       }
	  }
     } else {
	  if (SINGLE_PRECISION) {
	       if (p->sign == -1) {
		    CASE(p, fft, fftr4_);
	       } else {
		    CASE(p, fft, fftr4_un);
	       }
	  } else {
	       if (p->sign == -1) {
		    CASE(p, fft, fftr8_);
	       } else {
		    CASE(p, fft, fftr8_un);
	       }
	  }
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     void *in = p->in;
     void (*FFT)() = fft; /* cache */

     for (i = 0; i < iter; ++i) {
	  FFT(in);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
