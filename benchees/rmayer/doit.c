/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", NAME)
BENCH_DOC("package", "rmayer")
BENCH_DOC("author", "Ron Mayer")
BENCH_DOC("year", "1993")
BENCH_DOC("language", "C")
BENCH_DOC("email", "ramayer@geocities.com")
BENCH_DOC("url", "http://www.geocities.com/ResearchTriangle/8869/fft_summary.html")
BENCH_DOC("url-was-valid-on", "Wed Aug 15 22:39:07 EDT 2001")
BENCH_DOC("notes", "Complex data are stored in separate real/imag arrays.")
BENCH_DOC("notes", NOTE)
BENCH_DOC("notes", "Based on the Fast Hartley Transform (FHT).")
BENCH_DOC("notes", "The inverse transform is scaled.")
BENCH_DOC("copyright",
"Copyright 1988, 1993; Ron Mayer\n"
"\n"
"NOTE: This routine uses at least 2 patented algorithms, and may be\n"
"      under the restrictions of a bunch of different organizations.\n"
"      Although I wrote it completely myself; it is kind of a derivative\n"
"      of a routine I once authored and released under the GPL, so it\n"
"      may fall under the free software foundation's restrictions;\n"
"      it was worked on as a Stanford Univ project, so they claim\n"
"      some rights to it; it was further optimized at work here, so\n"
"      I think this company claims parts of it.  The patents are\n"
"      held by R. Bracewell (the FHT algorithm) and O. Buneman (the\n"
"      trig generator), both at Stanford Univ.\n"
"      If it were up to me, I'd say go do whatever you want with it;\n"
"      but it would be polite to give credit to the following people\n"
"      if you use this anywhere:\n"
"          Euler     - probable inventor of the fourier transform.\n"
"          Gauss     - probable inventor of the FFT.\n"
"          Hartley   - probable inventor of the hartley transform.\n"
"          Buneman   - for a really cool trig generator\n"
"          Mayer(me) - for authoring this particular version and\n"
"                      including all the optimizations in one package.\n"
"      Thanks,\n"
	  "      Ron Mayer; mayer@acuson.com\n")
END_BENCH_DOC

extern void fft(int n, double *r, double *i);
extern void ifft(int n, double *r, double *i);
extern void realfft(int n, double *r);
extern void realifft(int n, double *r);
extern void init_code(int n);

int can_do(struct problem *p)
{
     return (DOUBLE_PRECISION &&
	     p->rank == 1 &&
	     problem_power_of_two(p, 1) &&
	     p->n[0] > 2
	  );
}

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
#ifdef CALL_INIT_CODE
     init_code(p->n[0]);
#endif
}

void doit(int iter, struct problem *p)
{
     int n = p->size;
     int i;

     if (p->kind == PROBLEM_COMPLEX) {
	  double *re = (double *) p->in, *im = re + n;

	  if (p->sign < 0)
	       for (i = 0; i < iter; ++i)
		    fft(n, re, im);
	  else
	       for (i = 0; i < iter; ++i)
		    ifft(n, re, im);
     }
     else /* PROBLEM_REAL */ {
	  double *re = (double *) p->in;

	  if (p->sign < 0)
	       for (i = 0; i < iter; ++i)
		    realfft(n, re);
	  else
	       for (i = 0; i < iter; ++i)
		    realifft(n, re);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, 1);
}

void copy_c2c_from(struct problem *p, bench_complex *in)
{
     bench_real *x = (bench_real *) p->in;
     copy_c2ri(in, x, x + p->size, p->size);
}

void copy_c2c_to(struct problem *p, bench_complex *out)
{
     bench_real *x = (bench_real *) p->out;
     copy_ri2c(x, x + p->size, out, p->size);
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_1d_halfcomplex(p, out, 1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_halfcomplex(p, in, 1.0);
}

