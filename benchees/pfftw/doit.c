/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "pfftw")
BENCH_DOC("version", "0.03")
BENCH_DOC("author", "Matteo Frigo")
BENCH_DOC("email", "athena@fftw.org")
BENCH_DOC("year", "2001")
BENCH_DOC("language", "assembly")
BENCH_DOC("url", "http://www.fftw.org/download.html")
BENCH_DOC("url-was-valid-on", "Fri Mar 28 18:27:16 EST 2003")
BENCH_DOC("notes", "out of order")
BENCH_DOC("copyright", "pfftw is copyright 1999 by Matteo Frigo.  pfftw is distributed under the terms of the GNU General Public License.  See file COPYING for details.")
END_BENCH_DOC

#include "pfftw.h"

#if BENCHFFT_SINGLE
typedef fftw_complex_s C;
#else
typedef fftw_complex_d C;
#endif

static void (*fft)(C *);

int can_do(struct problem *p)
{
     return p->rank == 1 && problem_complex_power_of_two(p, 1) && 
	  p->n[0] <= 1024;
}

#define PASTEx(a,b) a ## b
#define PASTE(a,b) PASTEx(a,b)

#ifdef BENCHFFT_SINGLE
#  define PFFTW pfftw_s_
#  define PFFTWI pfftwi_s_
#else
#  define PFFTW pfftw_d_
#  define PFFTWI pfftwi_d_
#endif
#define PFFTWP(x) PASTE(PFFTW, x)

#define CASE(p, var, fft)			\
switch (p->n[0]) {				\
    case 2: var = PASTE(fft, 2); break;		\
    case 4: var = PASTE(fft, 4); break;		\
    case 8: var = PASTE(fft, 8); break;		\
    case 16: var = PASTE(fft, 16); break;      	\
    case 32: var = PASTE(fft, 32); break;      	\
    case 64: var = PASTE(fft, 64); break;      	\
    case 128: var = PASTE(fft, 128); break;    	\
    case 256: var = PASTE(fft, 256); break;    	\
    case 512: var = PASTE(fft, 512); break;    	\
    case 1024: var = PASTE(fft, 1024); break;	\
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

	  CASE(p, perm, PFFTWP(permutation_));

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

	  CASE(p, perm, PFFTWP(permutation_));

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
	  CASE(p, fft, PFFTW);
     } else {
	  CASE(p, fft, PFFTWI);
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     C *in = (C *) p->in;
     void (*FFT)(C*) = fft; /* cache */

     for (i = 0; i < iter; ++i) {
	  FFT(in);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
