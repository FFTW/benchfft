/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "mixfft")
BENCH_DOC("version", "0.5")
BENCH_DOC("author", "Jens Joergen Nielsen")
BENCH_DOC("year", "2000")
BENCH_DOC("language", "C")
BENCH_DOC("email", "jjn@get2net.dk")
BENCH_DOC("url", "http://hjem.get2net.dk/jjn/fft.htm")
BENCH_DOC("url-was-valid-on", "Mon Feb 10 01:19:57 EST 2003")
BENCH_DOC("notes", "Complex data are stored in separate real/imag arrays.")
BENCH_DOC("copyright",
"The source code in this packet is copyrighted material.\n"
"\n"
"Non-commercial use of the source code is free.\n"
"\n"
"A $100 fee must be paid if used commercially. Please contact me at\n"
"jjn@get2net.dk and register your copy. The $100 fee includes\n"
"up to one hour of assistance related to your use of the code. A\n"
"trial period of 14 days is allowed.\n"
"\n"
"If the code is used for professional (paid) research and development\n"
"for non-profit organisations like universities a reduced fee of $20\n"
"must be paid.\n"
"\n"
"The commercial license allows you to include the compiled code in a\n"
"product or to use the code on a regular basis. You are however NOT\n"
"allowed to sell the source code.\n"
"\n"
"Distribution of the complete unaltered package, including this\n"
"file, is free. This includes commercial CD's.\n")
END_BENCH_DOC

static const int maxPrimeFactor = 37;  /* must be same as in Mixfft.c */

int can_do(struct problem *p)
{
     return (DOUBLE_PRECISION &&
	     p->rank == 1 &&
	     p->kind == PROBLEM_COMPLEX &&
	     check_prime_factors(p->n[0], maxPrimeFactor) &&
	     !problem_in_place(p) &&
	     p->sign == -1
	  );
}

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
}

void problem_ccopy_from(struct problem *p, bench_complex *in)
{
     bench_real *x = (bench_real *) p->in;
     copy_c2ri(in, x, x + p->n[0], p->n[0]);
}

void problem_ccopy_to(struct problem *p, bench_complex *out)
{
     bench_real *x = (bench_real *) p->out;
     copy_ri2c(x, x + p->n[0], out, p->n[0]);
}

extern void fft(int n, double xRe[], double xIm[],
                double yRe[], double yIm[]);

void doit(int iter, struct problem *p)
{
     int n = p->n[0];
     bench_real *xr = (bench_real *) p->in;
     bench_real *xi = xr + n;
     bench_real *yr = (bench_real *) p->out;
     bench_real *yi = yr + n;
     int i;

     for (i = 0; i < iter; ++i) {
	  fft(n, xr, xi, yr, yi);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
