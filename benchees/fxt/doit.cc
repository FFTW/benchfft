/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

#include "fxt.h"

BEGIN_BENCH_DOC
BENCH_DOC("name", NAME)
BENCH_DOC("version", "26-August-2002")
BENCH_DOC("author", "J&ouml;rg Arndt")
BENCH_DOC("year", "2002")
BENCH_DOC("email", "arndt@jjj.de")
BENCH_DOC("url", "http://www.jjj.de/fxt/")
BENCH_DOC("url-was-valid-on", "Sun Jul 15 00:02:30 EDT 2001")
BENCH_DOC("copyright", 
"Copyright (c) 2001 Joerg Arndt\n"
"\n"
"This program is free software; you can redistribute it and/or modify\n"
"it under the terms of the GNU General Public License as published by\n"
"the Free Software Foundation; either version 2 of the License, or\n"
"(at your option) any later version.\n"
"\n"
"This program is distributed in the hope that it will be useful,\n"
"but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"GNU General Public License for more details.\n"
"\n"
"You should have received a copy of the GNU General Public License\n"
"along with this program; if not, write to the Free Software\n"
"Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA\n")
BENCH_DOC("language", "C++")
BENCH_DOC("notes", NOTES)
END_BENCH_DOC

#ifndef RANK_OK
#  define RANK_OK(r) ((r) == 1)
#endif

extern "C"
int can_do(struct problem *p)
{
     return (DOUBLE_PRECISION &&
	     p->kind == PROBLEM_COMPLEX &&
	     RANK_OK(p->rank) &&
	     problem_power_of_two(p, 1));
}

Complex *x = (Complex *) 0;
int m, m1, N, n0, n1;
ulong *md;

extern "C"
void setup(struct problem *p)
{
     int i, n = p->n[0];
#ifndef NO_Complex
     x = new Complex[n];
     for (i = 0; i < n; ++i) {
	  x[i] = Complex(0.0, 0.0);
     }
#endif
     md = new ulong[p->rank];
     N = 1;
     for (i = 0; i < p->rank; ++i) {
       md[i] = log_2(p->n[p->rank - 1 - i]);
       N *= p->n[i];
     }
     if (p->rank > 0) {
       n0 = p->n[0];
       m = md[p->rank - 1 - 0];
     }
     if (p->rank > 1) {
       n1 = p->n[1];
       m1 = md[p->rank - 1 - 1];
     }
}

#ifndef NO_Complex
extern "C"
void problem_ccopy_from(struct problem *p, bench_complex *in)
{
     int i;
     Complex *out = x;

     for (i = 0; i < N; ++i) {
	  out[i] = Complex(c_re(in[i]), c_im(in[i]));
     }
}

extern "C"
void problem_ccopy_to(struct problem *p, bench_complex *out)
{
     int i;
     Complex *in = x;

     for (i = 0; i < N; ++i) {
	  c_re(out[i]) = real(in[i]);
	  c_im(out[i]) = imag(in[i]);
     }
}
#else /* NO_Complex */
void problem_ccopy_from(struct problem *p, bench_complex *in)
{
     bench_real *x = (bench_real *) p->in;
     copy_c2ri(in, x, x + N, N);
}

void problem_ccopy_to(struct problem *p, bench_complex *out)
{
     bench_real *x = (bench_real *) p->in;
     copy_ri2c(x, x + N, out, N);
}
#endif /* NO_Complex */

extern "C"
void doit(int iter, struct problem *p)
{
     int i;
     int is = p->sign;
#ifdef NO_Complex
     int rank= p->rank;
     double *inr = (double *) p->in;
     double *ini = inr + N;
#endif

     for (i = 0; i < iter; ++i) {
	  DOIT_FFT;
     }
}

extern "C"
void done(struct problem *p)
{
     UNUSED(p);
     delete md;
#ifndef NO_Complex
     delete x;
#endif
}
