/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
#ifdef BLOODWORTH_FHT
  BENCH_DOC("name", "bloodworth-fht")
#else
  BENCH_DOC("name", "bloodworth")
#endif
BENCH_DOC("package", "bloodworth")
BENCH_DOC("author", "Carey E. Bloodworth")
BENCH_DOC("year", "1998")
BENCH_DOC("language", "C")
BENCH_DOC("email", "cbloodworth@juno.com")
BENCH_DOC("notes", "Received in personal communication from the author.")
BENCH_DOC("notes", "Slightly modified for benchFFT 2.0 (some subroutine names changed, and code slightly altered to work with both single and double precision...original code was double precision, I think).")
#ifdef BLOODWORTH_FHT
  BENCH_DOC("notes", "Uses a Fast Hartley Transform (FHT) algorithm.")
#else
  BENCH_DOC("notes", "The backwards transform is scaled.")
#endif
BENCH_DOC("copyright",
"The code is original to me.  Although I have copyrighted it, you can "
"freely use it for any reason provided my name remains.  The only reason "
"it's not PD is because I put quite a bit of work into it and I just like "
"the idea of my name being in code that somebody else is using.")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (p->rank == 1 &&
	     problem_power_of_two(p, 1)
#ifdef BLOODWORTH_FHT
	     && p->n[0] >= 8
	     && p->kind == PROBLEM_REAL
#endif
	  );
}

extern void Bloodworth_Q2_InitFFT(unsigned int n);
extern void Bloodworth_Q2_FwdFFT(bench_complex *data, int n);
extern void Bloodworth_Q2_RevFFT(bench_complex *data, int n);
extern void Bloodworth_Q2_FwdRealFFT(bench_real *data, int n);
extern void Bloodworth_Q2_RevRealFFT(bench_real *data, int n);

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
     Bloodworth_Q2_InitFFT(p->size);
}

void doit(int iter, struct problem *p)
{
     int n = p->size;
     int i;

#ifndef BLOODWORTH_FHT
     if (p->kind == PROBLEM_COMPLEX) {
	  bench_complex *in = (bench_complex *) p->in;

	  if (p->sign > 0)
	       for (i = 0; i < iter; ++i)
		    Bloodworth_Q2_FwdFFT(in, n);
	  else
	       for (i = 0; i < iter; ++i)
		    Bloodworth_Q2_RevFFT(in, n);
     }
     else /* PROBLEM_REAL */
#endif
     {
	  bench_real *in = (bench_real *) p->in;

	  if (p->sign < 0)
	       for (i = 0; i < iter; ++i)
		    Bloodworth_Q2_FwdRealFFT(in, n);
	  else
	       for (i = 0; i < iter; ++i)
		    Bloodworth_Q2_RevRealFFT(in, n);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
#ifndef BLOODWORTH_FHT
     unnormalize(p, out, p->kind == PROBLEM_COMPLEX ? -1 : 1);
#endif
}

void copy_h2c(struct problem *p, bench_complex *out)
{
#ifdef BLOODWORTH_FHT
     copy_h2c_1d_halfcomplex(p, out, 1.0);
#else
     copy_h2c_1d_packed(p, out, 1.0);
#endif
}

void copy_c2h(struct problem *p, bench_complex *in)
{
#ifdef BLOODWORTH_FHT
     copy_c2h_1d_halfcomplex(p, in, 1.0);
#else
     copy_c2h_1d_packed(p, in, 1.0);
#endif
}

