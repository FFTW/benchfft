/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "sorensen-ctfftsr")
BENCH_DOC("author", "Henrik Sorensen")
BENCH_DOC("year", "1984")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("notes", "Complex data are stored in separate real/imag arrays.")
BENCH_DOC("bibitem", 
	  "Sorensen, Heideman, Burrus: On computing the split-radix FFT,\n"
	  "IEEE Tran. ASSP, Vol. ASSP-34, No. 1, pp. 152-156 Feb. 1986")
BENCH_DOC("notes",
	  "no inverse is provided, although one can easily swap the real and"
	  " imaginary pointers.")
END_BENCH_DOC

#define _CTFFTSR F77_FUNC(ctfftsr, CTFFTSR)
#define _TINIT F77_FUNC(tinit, TINIT)

extern void _CTFFTSR();
extern void _TINIT();

int can_do(struct problem *p)
{
     return (SINGLE_PRECISION && p->rank == 1 && p->sign == -1 &&
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
static bench_real *CT1, *CT3, *ST1, *ST3;
static int *ITAB;

void setup(struct problem *p)
{
     int n;
     BENCH_ASSERT(can_do(p));
     
     n = p->n[0];
     M = log_2(n);
     CT1 = bench_malloc(n / 8 * sizeof(bench_real));
     CT3 = bench_malloc(n / 8 * sizeof(bench_real));
     ST1 = bench_malloc(n / 8 * sizeof(bench_real));
     ST3 = bench_malloc(n / 8 * sizeof(bench_real));
     ITAB = bench_malloc((1 << (M / 2 + 1)) * sizeof(bench_real));
     _TINIT(&M, CT1, CT3, ST1, ST3, ITAB);
}

void doit(int iter, struct problem *p)
{
     int i, n = p->n[0], m = M;
     bench_real *x = p->in;
     bench_real *y = x + n;
     bench_real *ct1 = CT1, *st1 = ST1;
     bench_real *ct3 = CT3, *st3 = ST3;
     int *itab = ITAB;

     for (i = 0; i < iter; ++i) {
	  _CTFFTSR(x, y, &m, ct1, ct3, st1, st3, itab);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
