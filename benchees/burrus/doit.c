/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "burrus-sffteu")
BENCH_DOC("author", "C. S. Burrus")
BENCH_DOC("year", "1984")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("notes", "Complex data are stored in separate real/imag arrays.")
BENCH_DOC("bibitem", "Electronics Letters, January 5, 1984")
BENCH_DOC("url", 
	  "http://ourworld.compuserve.com/homepages/steve_kifowit/sffteu.txt")
BENCH_DOC("url-was-valid-on", "Mon Jul 16 00:03:31 EDT 2001")
END_BENCH_DOC

#define _SFFTEU F77_FUNC(sffteu, SFFTEU)

extern void _SFFTEU();

int can_do(struct problem *p)
{
     return (SINGLE_PRECISION && p->rank == 1 &&
	     problem_complex_power_of_two(p, 1));
}


void problem_ccopy_from(struct problem *p, bench_complex *in)
{
     BENCH_ASSERT(p->kind == PROBLEM_COMPLEX && p->rank == 1);
     copy_c2ri(in, p->in, ((bench_real *)p->in) + p->n[0],
	       p->n[0]);
}

void problem_ccopy_to(struct problem *p, bench_complex *out)
{
     BENCH_ASSERT(p->kind == PROBLEM_COMPLEX && p->rank == 1);
     copy_ri2c(p->out, ((bench_real *)p->out) + p->n[0], out,
	       p->n[0]);

     unnormalize(p, out, 1);
}

static int M;
void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
     M = log_2(p->n[0]);
}

void doit(int iter, struct problem *p)
{
     int i, n = p->n[0], m = M, itype = -p->sign;
     bench_real *x = p->in;
     bench_real *y = x + n;

     for (i = 0; i < iter; ++i) {
	  _SFFTEU(x, y, &n, &m, &itype);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
