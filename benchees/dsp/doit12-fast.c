/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "dsp79-FAST")
BENCH_DOC("package", "Programs for Digital Signal Processing")
BENCH_DOC("author", "G. D. Bergland and M. T. Dolan")
BENCH_DOC("year", "1979")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("notes", "Subroutine FAST from DSP Section 1.2")
BENCH_DOC("notes", "The backward transform is scaled")
BENCH_DOC("bibitem", "IEEE DSP Committee, Programs for Digital Signal Processing (IEEE Press, 1979). ISBN 0-471-05962-5. (Out of Print)")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (SINGLE_PRECISION &&
	     p->rank == 1 &&
	     problem_real_power_of_two(p, 1) &&
	     p->n[0] <= 32768);
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_unpacked(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_unpacked(p, in, -1.0);
}

void copy_r2c(struct problem *p, bench_complex *out)
{
     unsigned int i;
     unsigned int n = p->size;
     bench_real *pout = p->out;
     bench_real dn = (bench_real)n;

     for (i = 0; i < n; ++i) {
	  c_re(out[i]) = pout[i] * dn;
	  c_im(out[i]) = 0.0;
     }
}

extern void F77_FUNC(fast, FAST)();
extern void F77_FUNC(fsst, FSST)();

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
}

void doit(int iter, struct problem *p)
{
     int n = p->n[0];
     bench_real *in = p->in;
     int i;

     if (p->sign == -1) {
	  for (i = 0; i < iter; ++i) {
	       F77_FUNC(fast, FAST)(in, &n);
	  }
     } else {
	  for (i = 0; i < iter; ++i) {
	       F77_FUNC(fsst, FSST)(in, &n);
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
