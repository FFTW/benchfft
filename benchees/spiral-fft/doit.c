/* this program is in the public domain */

#include <stdio.h>
#include "bench-user.h"
#undef c_re
#undef c_im
#include "spiral_fft.h"

BEGIN_BENCH_DOC
BENCH_DOC("name", "spiral-egner-fft")
BENCH_DOC("version", "1")
BENCH_DOC("author", "Gavin Haentjens")
BENCH_DOC("email", "gh2r@andrew.cmu.edu")
BENCH_DOC("author", "Adrian Sox")
BENCH_DOC("email", "asox@andrew.cmu.edu")
BENCH_DOC("year", "2000")
BENCH_DOC("language", "C")
BENCH_DOC("url", "http://www.ece.cmu.edu/~spiral/fft.html")
BENCH_DOC("url-was-valid-on", "Sun Sep  1 19:30:13 EDT 2002")
BENCH_DOC("notes", "Design is based on an FFT package by Sebastian Egner.")
BENCH_DOC("notes", "Backwards transform is forward transform plus two extra passes to scale and reverse the output.")
BENCH_DOC("notes", "Uses codelets (hard-coded FFTs of small sizes, with or without premultiplying twiddle factors) from FFTW 2.1.1.")
BENCH_DOC("copyright",
"Copyright (c) 2000 Carnegie Mellon University\n"
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
END_BENCH_DOC

/* max FFT is 2^(MAXSIZE-1); MAXSIZE is set in spiral_fft/configure.in */
#define MAXSIZE 20

int can_do(struct problem *p)
{
     return (DOUBLE_PRECISION &&
	     p->rank == 1 &&
	     p->kind == PROBLEM_COMPLEX &&
	     p->n[0] < (1 << MAXSIZE) &&
	     problem_power_of_two(p, 1));
}

int m = 0;
fft_t *tree = 0;

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
     m = log_2(p->n[0]);
     BENCH_ASSERT(tree = get_tree(m));
     if (verbose > 1) {
	  char *s;
	  s = fft_print(2, tree);
	  printf("%s\n", s);
     }
}

#define N 1

void doit(int iter, struct problem *p)
{
     fft_value *x = (fft_value *) p->in;
     int i;

     if (p->sign < 0)
	  for (i = 0; i < iter; ++i) {
	       fft_apply(tree, 1, x);
	  }
     else
	  for (i = 0; i < iter; ++i) {
	       fft_apply_inverse(tree, 1, x);
	  }
}

void done(struct problem *p)
{
     UNUSED(p);
}
