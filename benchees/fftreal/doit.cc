/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "fftreal")
BENCH_DOC("package", "FFTReal")
BENCH_DOC("year", "2001")
BENCH_DOC("version", "1.02")
BENCH_DOC("author", "Laurent de Soras")
BENCH_DOC("email", "ldesoras@club-internet.fr")
BENCH_DOC("copyright",
	  "Source code may be freely used for any purpose, including commercial applications. Programs must display in their 'About' dialog-box (or documentation) a text telling they use these routines by Laurent de Soras. Modified source code can be distributed, but modifications must be clearly indicated.")
BENCH_DOC("language", "C++")
BENCH_DOC("notes", "Code is in single precision, but double precision version is created as documented by the author, by changing the flt_t typedef in the FFTReal.h header file")
BENCH_DOC("notes", "According to the README: There is no official web site where to get these files. If you have them, that's all to the good for you.")
BENCH_DOC("url", "http://www.musicdsp.org/FFTReal.zip")
BENCH_DOC("url-was-valid-on", "Tue Apr 15 17:01:59 EDT 2003")
END_BENCH_DOC

#include "FFTReal.h"

int can_do(struct problem *p)
{
  return (problem_power_of_two(p, 0) &&
	  p->kind == PROBLEM_REAL
	  && p->rank == 1);
}

void copy_h2c(struct problem *p, bench_complex *out)
{
  bench_real *h = (bench_real *) p->out;
  int i;
  int n = p->n[0];

  c_re(out[0]) = h[0];
  for (i = 1; i <= n/2; ++i)
    c_re(out[n-i]) = c_re(out[i]) = h[i];
  c_im(out[0]) = c_im(out[n/2]) = 0.0;
  for (i = 1; i < n/2; ++i)
    c_im(out[n-i]) = -(c_im(out[i]) = -h[n/2 + i]);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
  bench_real *h = (bench_real *) p->in;
  int i;
  int n = p->n[0];

  for (i = 0; i <= n/2; ++i)
    h[i] = c_re(in[i]);
  for (i = 1; i < n/2; ++i)
    h[n/2 + i] = -c_im(in[i]);
}

FFTReal *fft = (FFTReal *) 0;

void setup(struct problem *p)
{
  fft = new FFTReal(p->n[0]);
}

void doit(int iter, struct problem *p)
{
     int i;
     const bench_real *in = (const bench_real *) p->in;
     bench_real *out = (bench_real *) p->out;

     if (p->sign == -1) {
       for (i = 0; i < iter; ++i) {
	 fft->do_fft(out, in);
       }
     } else {
       for (i = 0; i < iter; ++i) {
	 fft->do_ifft(in, out);
       }
     }
}

void done(struct problem *p)
{
  delete fft;
}
