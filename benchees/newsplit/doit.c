/* this program is in the public domain */

#include "bench-user.h"
#include "split.h"

BEGIN_BENCH_DOC
BENCH_DOC("name", NAME)
BENCH_DOC("author", "Steven G. Johnson")
BENCH_DOC("year", "2005")
BENCH_DOC("language", "C")
BENCH_DOC("email", "stevenj@alum.mit.edu")
BENCH_DOC("notes", NOTE)
BENCH_DOC("notes", "This is a toy implementation to demonstrate a new version of the split-radix FFT with reduced arithmetic complexity.  It is not intended for actual use (being very slow, and lacking an inverse transform).")
BENCH_DOC("copyright",
"Copyright (c) 2005 Massachusetts Institute of Technology\n"
"\n"
"Permission is hereby granted, free of charge, to any person obtaining\n"
"a copy of this software and associated documentation files (the\n"
"\"Software\"), to deal in the Software without restriction, including\n"
"without limitation the rights to use, copy, modify, merge, publish,\n"
"distribute, sublicense, and/or sell copies of the Software, and to\n"
"permit persons to whom the Software is furnished to do so, subject to\n"
"the following conditions:\n"
"\n"
"The above copyright notice and this permission notice shall be\n"
"included in all copies or substantial portions of the Software.\n"
"\n"
"THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND,\n"
"EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF\n"
"MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.\n"
"IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY\n"
"CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,\n"
"TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE\n"
"SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (p->rank == 1 &&
	     p->sign == -1 &&
	     p->kind == PROBLEM_COMPLEX &&
	     problem_power_of_two(p, 0));
}

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));

     if (verbose >= 2) {
	  adds = muls = twids = 0;
	  doit(1, p);
	  printf("adds = %d, muls = %d, flops = %d\n", adds, muls, adds+muls);
	  printf("real twiddle loads = %d\n", twids);
     }
#ifndef OLDSPLIT
     if (verbose >= 2) {
	  int k, n = p->n[0];
	  int kmin = 0;
	  double smin = 1.0;
	  for (k = 0; k < n/4; ++k) {
	       double sk = scale(n, k);
	       if (sk < smin && fabs(sk - smin) > 1e-12 * fabs(smin)) {
		    kmin = k;
		    smin = sk;
	       }
	  }
	  printf("minimum scale(%d, %d) = %g\n", n, kmin, smin);
     }
#endif
}

void doit(int iter, struct problem *p)
{
     int n = p->n[0];
     C *in = (C *) p->in;
     C *out = (C *) p->out;
     int i;

     for (i = 0; i < iter; ++i) {
	  FFT(n, in, in+1, 1, out, 1);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
