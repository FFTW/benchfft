/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "krukar")
BENCH_DOC("author", "Richard H. Krukar")
BENCH_DOC("year", "1990")
BENCH_DOC("language", "C")
BENCH_DOC("email", "krukar@ectopia.com")
BENCH_DOC("url", "http://risc2.numis.nwu.edu/ftp/pub/transforms/ffts_in_C.tar.gz")
BENCH_DOC("url-was-valid-on", "Sat Aug 31 16:36:36 EDT 2002")
BENCH_DOC("notes", "The backward transform is scaled")
BENCH_DOC("copyright",
"copyright 1990  Richard H. Krukar all rights reserved\n"
"\n"
"Permission granted to buy, sell, or steal this software is granted.\n"
"The author retains the right to distribute this software freely, but\n"
"is not responsible for it's quality or maintainance.")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (SINGLE_PRECISION &&
	     p->rank == 1 &&
	     p->kind == PROBLEM_COMPLEX &&
	     problem_power_of_two(p, 1) &&
	     p->n[0] <= 4096);
}

typedef void (*krukar_fft_t)(bench_complex *x, int n);
extern void fft(bench_complex *x, int n);
extern void ifft(bench_complex *x, int n);

krukar_fft_t krukar_fft;

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
     if (p->sign == +1)
	  krukar_fft = ifft;
     else
	  krukar_fft = fft;
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, 1);
}

void doit(int iter, struct problem *p)
{
     int n = p->n[0];
     bench_complex *x = (bench_complex *) p->in;
     int i;

     for (i = 0; i < iter; ++i) {
	  krukar_fft(x, n);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
