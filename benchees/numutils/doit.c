/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "numutils")
BENCH_DOC("author", "Carl Hein")
BENCH_DOC("year", "2002")
BENCH_DOC("language", "C")
BENCH_DOC("email", "chein@atl.lmco.com")
BENCH_DOC("url", "http://www.atl.external.lmco.com/proj/csim/xgraph/numutil/")
BENCH_DOC("url-was-valid-on", "Sun Dec  1 21:06:34 EST 2002")
BENCH_DOC("notes", "The forward transform is scaled")
BENCH_DOC("copyright",
"Copyright (C) 2002 CHein\n"
"\n"
"This program is free software; you can redistribute it and/or modify\n"
"it under the terms of the GNU General Public License as published by\n"
"the Free Software Foundation, Version 2 of the GPL.  This program is\n"
"distributed in the hope that it will be useful, but WITHOUT ANY\n"
"WARRANTY; without even the implied warranty of MERCHANTABILITY or\n"
"FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License\n"
"for more details.  You should have received a copy of the GNU General\n"
"Public License along with this program; if not, write to the Free\n"
"Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA\n"
"02111-1307.")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (DOUBLE_PRECISION &&
	     p->rank == 1 &&
	     p->kind == PROBLEM_COMPLEX &&
	     !problem_in_place(p));
}

extern void FFT(bench_complex *x, bench_complex *Y, int N, int direction, 
		bench_real *SN, bench_real *CN);

bench_real *SN = 0, *CN = 0; /* trig. tables */

void setup(struct problem *p)
{
     int n = p->n[0], dir = -p->sign;
     bench_complex *x = (bench_complex *) p->in;
     bench_complex *y = (bench_complex *) p->out;
     int i;

     BENCH_ASSERT(can_do(p));
     SN = (bench_real *) bench_malloc(sizeof(bench_real) * n);
     CN = (bench_real *) bench_malloc(sizeof(bench_real) * n);
     FFT(x, y, n, dir, SN, CN); /* initialize trig. tables on first call */
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, -1);
}

void doit(int iter, struct problem *p)
{
     int n = p->n[0], dir = -p->sign;
     bench_complex *x = (bench_complex *) p->in;
     bench_complex *y = (bench_complex *) p->out;
     int i;

     for (i = 0; i < iter; ++i) {
	  FFT(x, y, n, dir, SN, CN);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(CN);
     bench_free(SN);
}
