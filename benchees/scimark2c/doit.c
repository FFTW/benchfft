/* this program is in the public domain */

#include "bench-user.h"

#include "FFT.h"

BEGIN_BENCH_DOC

BENCH_DOC("name", "scimark2c")
BENCH_DOC("package", "SciMark")
BENCH_DOC("version", "2.0")
BENCH_DOC("language", "C")
BENCH_DOC("year", "2000")
BENCH_DOC("author", "Roldan Pozo")
BENCH_DOC("email", "pozo@nist.gov")
BENCH_DOC("author", "Bruce Miller")
BENCH_DOC("email", "bruce.miller@nist.gov")
BENCH_DOC("url", "http://math.nist.gov/scimark2/")
BENCH_DOC("url-was-valid-on", "Tue May  3 22:33:43 EDT 2005")
BENCH_DOC("notes", "C version of SciMark Java benchmark.")
BENCH_DOC("notes", "The backwards transform is scaled.")

END_BENCH_DOC

int can_do(struct problem *p)
{
     return (DOUBLE_PRECISION &&
	     p->rank == 1 &&
	     p->kind == PROBLEM_COMPLEX &&
	     problem_power_of_two(p, 1)
	  );
}

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
}

void done(struct problem *p)
{
     UNUSED(p);
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, +1);
}

void doit(int iter, struct problem *p)
{
     double *data = (double *) p->in;
     int N = p->size * 2;
     int i;

     if (p->sign < 0)
	  for (i = 0; i < iter; ++i)
	       FFT_transform(N, data);
     else
	  for (i = 0; i < iter; ++i)
	       FFT_inverse(N, data);
}
