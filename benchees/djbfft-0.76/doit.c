/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "djbfft-0.76")
BENCH_DOC("author", "D. J. Bernstein")
BENCH_DOC("year", "1999")
BENCH_DOC("email", "djb@pobox.com")
BENCH_DOC("url", "http://pobox.com/~djb/djbfft.html")
BENCH_DOC("language", "C")
BENCH_DOC("notes",
	  "djbfft considers sign = +1 to be the forward transform")
BENCH_DOC("conf-opt", CONF_OPT)
BENCH_DOC("auto-opt", AUTO_OPT)
END_BENCH_DOC

#include "fftc4.h"
#include "fftc8.h"
#include "fftfreq.h"

static void (*fft)();

int can_do(struct problem *p)
{
     return p->rank == 1 &&
	  problem_complex_power_of_two(p, 1) && 
	  p->n[0] <= 8192;
}

void problem_ccopy_from(struct problem *p, bench_complex *in)
{
     if (p->sign == 1)
	  /* in-order input for backward transforms */
	  cacopy(in, p->in, p->n[0]);
     else {
	  unsigned int i;
	  bench_complex *pin = p->in;
	  for (i = 0; i < p->n[0]; ++i) {
	       pin[i] = in[fftfreq_c(i, p->n[0])];
	  }
     }
}

void problem_ccopy_to(struct problem *p, bench_complex *out)
{
     if (p->sign == 1) {
	  unsigned int i;
	  bench_complex *pout = p->out;
	  for (i = 0; i < p->n[0]; ++i) {
	       out[fftfreq_c(i, p->n[0])] = pout[i];
	  }
     } else {
	  /* in-order output for forward transforms */
	  cacopy(p->out, out, p->n[0]);
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
     if (sizeof(bench_real) == sizeof(float)) {
	  if (p->sign == 1) {
	       CASE(p, fft, fftc4_);
	  } else {
	       CASE(p, fft, fftc4_un);
	  }
     } else {
	  if (p->sign == 1) {
	       CASE(p, fft, fftc8_);
	  } else {
	       CASE(p, fft, fftc8_un);
	  }
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     bench_complex *in = p->in;
     void (*FFT)() = fft; /* cache */

     for (i = 0; i < iter; ++i) {
	  FFT(in);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
