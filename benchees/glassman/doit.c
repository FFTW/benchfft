/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "glassman")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("url", "http://www.jjj.de/fft/glassman-fft.f")
BENCH_DOC("url-was-valid-on", "Sun Nov 17 22:40:20 EST 2002")
BENCH_DOC("bibitem",
	  "W. E. Ferguson, Jr., A simple derivation of Glassman general-n fast Fourier transform, Comput. and Math. with Appls. 8 (6), 401-411 (1982). Also in Report AD-A083 811, NTIS, Dec. 1979.")
BENCH_DOC("bibitem", 
	  "J. A. Glassman, A generalization of the fast Fourier transform, IEEE Trans. on Computers 19 (2), 105-116 (1970).")
BENCH_DOC("notes", "The backward transform is scaled.")
BENCH_DOC("notes", "Uses an unstable trigonometric iteration.")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (SINGLE_PRECISION &&
	     p->rank == 1 &&
	     p->kind == PROBLEM_COMPLEX &&
	     problem_in_place(p)
	  );
}

#define SPCFFT F77_FUNC(spcfft,SPCFFT)

extern void SPCFFT(bench_complex *x, int *n, int *isign, bench_complex *work, float *interp);

static float interp = 1.0;
bench_complex *work;

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
     work = (bench_complex *) bench_malloc(sizeof(bench_complex) * p->size);
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, 1.0);
}

void doit(int iter, struct problem *p)
{
     int i;
     bench_complex *x = (bench_complex *) p->in;
     int isign = p->sign, n = p->size;

     for (i = 0; i < iter; ++i) {
	  SPCFFT(x, &n, &isign, work, &interp);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(work);
}

