/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "dsp79-wfta")
BENCH_DOC("package", "Programs for Digital Signal Processing")
BENCH_DOC("author", "James H. McClellan and Hamid Nawab")
BENCH_DOC("year", "1979")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("notes", "Complex data are stored in separate real/imag arrays.")
BENCH_DOC("notes", "The backward transform is scaled")
BENCH_DOC("notes", "DSP Section 1.7")
BENCH_DOC("notes", "Modified by Matteo Frigo to remove assumption that local variables are SAVEd")
BENCH_DOC("bibitem", "IEEE DSP Committee, Programs for Digital Signal Processing (IEEE Press, 1979). ISBN 0-471-05962-5. (Out of Print)")
END_BENCH_DOC

int can_do(struct problem *p)
{
     unsigned int n;

     if (!(SINGLE_PRECISION && p->rank == 1))
	  return 0;
     if (!problem_in_place(p))
	  return 0;
     if (p->kind != PROBLEM_COMPLEX)
	  return 0;

     n = p->n[0];

     if (n == 0 || n > 5040)
	  return 0;

     if (n % 16 == 0) n /= 16;
     else if (n % 8 == 0) n /= 8;
     else if (n % 4 == 0) n /= 4;
     else if (n % 2 == 0) n /= 2;

     if (n % 9 == 0) n /= 9;
     else if (n % 3 == 0) n /= 3;

     if (n % 7 == 0) n /= 7;

     if (n % 5 == 0) n /= 5;

     return (n == 1);
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

extern void F77_FUNC(wfta, WFTA)();

void setup(struct problem *p)
{
     int n = p->n[0];
     int invrs = (p->sign == 1);
     bench_real *rin = p->in;
     bench_real *iin = rin + n;
     int init = 0;
     int ierr = 0;

     BENCH_ASSERT(can_do(p));
     F77_FUNC(wfta, WFTA)(rin, iin, &n, &invrs, &init, &ierr);
}

void doit(int iter, struct problem *p)
{
     int n = p->n[0];
     int invrs = (p->sign == 1);
     bench_real *rin = p->in;
     bench_real *iin = rin + n;
     int init = 1;
     int ierr = 0;
     int i;

     for (i = 0; i < iter; ++i) {
	  F77_FUNC(wfta, WFTA)(rin, iin, &n, &invrs, &init, &ierr);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
