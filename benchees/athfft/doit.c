/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "athfft")
BENCH_DOC("author", "Matteo Frigo")
BENCH_DOC("email", "athena@fftw.org")
BENCH_DOC("year", "2001")
BENCH_DOC("language", "C")
BENCH_DOC("notes", "out of order")
END_BENCH_DOC

#include "athfft.h"

static void (*fft)();

int can_do(struct problem *p)
{
     return p->rank == 1 && problem_complex_power_of_two(p, 1) && 
	  p->n[0] <= 1024;
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
}

void copy_c2c_from(struct problem *p, bench_complex *in)
{
     if (p->sign == -1)
	  /* in-order input for forward transforms */
	  cacopy(in, p->in, p->n[0]);
     else {
	  unsigned int i;
	  bench_complex *pin = p->in;
	  int (*perm)(int);

	  CASE(p, perm, ath_permutation_);

	  for (i = 0; i < p->n[0]; ++i) {
	       pin[perm(i)] = in[i];
	  }
     }
}

void copy_c2c_to(struct problem *p, bench_complex *out)
{
     if (p->sign == -1) {
	  unsigned int i;
	  bench_complex *pout = p->out;
	  int (*perm)(int);

	  CASE(p, perm, ath_permutation_);

	  for (i = 0; i < p->n[0]; ++i) {
	       out[i] = pout[perm(i)];
	  }
     } else {
	  /* in-order output for backward transforms */
	  cacopy(p->out, out, p->n[0]);
     }
}


void setup(struct problem *p)
{
     if (p->sign == -1) {
	  CASE(p, fft, athfft_);
     } else {
	  CASE(p, fft, athffti_);
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     void *in = p->in;
     void (*FFT)() = fft; /* cache */

     for (i = 0; i < iter; ++i) {
	  FFT(in, 1, 0);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
