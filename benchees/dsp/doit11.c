/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "dsp79-rader")
BENCH_DOC("author", "C. M. Rader, after Brenner")
BENCH_DOC("year", "1979, after Brenner 1967")
BENCH_DOC("language", "FORTRAN")
BENCH_DOC("notes", "The backward transform is scaled")
BENCH_DOC("notes", "DSP Section 1.1")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (SINGLE_PRECISION &&
	     p->rank == 1 &&
	     problem_complex_power_of_two(p, 1) &&
	     p->n[0] <= 32768);
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, 1);
}

extern void F77_FUNC(fourea, FOUREA)();

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
}

void doit(int iter, struct problem *p)
{
     int n = p->n[0];
     int sign = p->sign;
     bench_complex *in = p->in;
     int i;

     for (i = 0; i < iter; ++i) {
	  F77_FUNC(fourea, FOUREA)(in, &n, &sign);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
