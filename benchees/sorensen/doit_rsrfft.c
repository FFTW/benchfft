/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "sorensen-rsrfft")
BENCH_DOC("author", "Henrik Sorensen")
BENCH_DOC("year", "1985")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("bibitem", 
	  "Sorensen, Jones, Heideman, Burrus: Real-valued fast\n" 
	  "Fourier transform algorithms, IEEE Tran. ASSP,\n"
	  "Vol. ASSP-35, No. 6, pp. 849-864, June 1987")
END_BENCH_DOC

#define _RSRFFT F77_FUNC(rsrfft, RSRFFT)
#define _IRSRFFT F77_FUNC(irsrfft, IRSRFFT)

extern void _RSRFFT();
extern void _IRSRFFT();

int can_do(struct problem *p)
{
     return (SINGLE_PRECISION && p->rank == 1 &&
	     p->n[0] <= 16384 &&  /* becomes unstable for large N */
#ifdef INVERSE
	     p->sign == 1
#else
	     p->sign == -1
#endif
	     && problem_real_power_of_two(p, 1));
}


void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, 1);
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_1d_halfcomplex(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_halfcomplex(p, in, -1.0);
}

static int M;
void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
     M = log_2(p->n[0]);
}

void doit(int iter, struct problem *p)
{
     int i, m = M;
     bench_real *x = p->in;

     for (i = 0; i < iter; ++i) {
	  /* cannot link both RSRFFT and IRSRFFT in the same program
	     because of redefinition of RBITREV */
#ifdef INVERSE
	  _IRSRFFT(x, &m);
#else
	  _RSRFFT(x, &m);
#endif
     }

}

void done(struct problem *p)
{
     UNUSED(p);
}
