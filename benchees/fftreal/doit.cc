/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "fftreal")
BENCH_DOC("package", "FFTReal")
BENCH_DOC("year", "2005")
BENCH_DOC("version", "2.00")
BENCH_DOC("author", "Laurent de Soras")
BENCH_DOC("email", "ldesoras@club-internet.fr")
BENCH_DOC("copyright",
	  "(c) Laurent de Soras <laurent.de.soras@club-internet.fr>\n"
"Object Pascal port (c) Frederic Vanmol <frederic@fruityloops.com>\n"
"\n"
"This library is free software; you can redistribute it and/or\n"
"modify it under the terms of the GNU Lesser General Public\n"
"License as published by the Free Software Foundation; either\n"
"version 2.1 of the License, or (at your option) any later version.\n"
"\n"
"This library is distributed in the hope that it will be useful,\n"
"but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU\n"
"Lesser General Public License for more details.\n"
"\n"
"You should have received a copy of the GNU Lesser General Public\n"
"License along with this library; if not, write to the Free Software\n"
"Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA")
BENCH_DOC("language", "C++")
BENCH_DOC("url", "http://ldesoras.free.fr/prod.html")
BENCH_DOC("url-was-valid-on", "Sat Oct 25 14:40:21 EDT 2008")
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

FFTReal<bench_real> *fft = NULL;

void setup(struct problem *p)
{
  fft = new FFTReal<bench_real>(p->n[0]);
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
