/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "harm")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("url", "http://risc2.numis.nwu.edu/ftp/pub/transforms/")
BENCH_DOC("url-was-valid-on", "Sat Aug 31 16:37:44 EDT 2002")
BENCH_DOC("bibitem", 
	  "May be related to radix-4 PK HARM subroutine mentioned in "
	  "J. W. Cooley and P. A. W. Lewis and P. D. Welch, "
	  "The Fast Fourier Transform Algorithm and Its Applications "
	  "(IBM Research, 1967).")

BENCH_DOC("copyright", "THIS ROUTINE IS PUBLIC DOMAIN")
BENCH_DOC("notes", "The backward transform is scaled")
END_BENCH_DOC

int ok_sizes(int rank, unsigned int *n)
{
     int i;
     for (i = 0; i < rank; ++i)
	  if ((1 << 17) < n[i] || 2 >= n[i]) /* suppress buggy sizes */
	       return 0;
     return 1;
}

int can_do(struct problem *p)
{
     return (p->rank >= 1 && p->rank <= 3 &&
	     p->kind == PROBLEM_COMPLEX &&
	     power_of_two(p->size) &&
	     problem_in_place(p) &&
	     ok_sizes(p->rank, p->n)
	  );
}

#ifdef BENCHFFT_SINGLE
#define HARM_F77 F77_FUNC(harm,HARM)
#else
#define HARM_F77 F77_FUNC(harmd,HARMD)
#endif

extern void HARM_F77(bench_complex *, int *, int *, bench_real *,
		     int *, int *);

int m[3], *inv;
bench_real *s;

void setup(struct problem *p)
{
     int ifset = 0, iferr = 0;
     unsigned int i, maxdim = 0;
     bench_complex *a = (bench_complex *) p->in;
     BENCH_ASSERT(can_do(p));

     m[0] = m[1] = m[2] = 0;
     for (i = 0; i < p->rank; ++i) {
	  m[i] = log_2(p->n[p->rank - 1 - i]);
	  if (p->n[i] > maxdim)
	       maxdim = p->n[i];
     }
     inv = (int*) bench_malloc(maxdim * sizeof(int));
     s = (bench_real*) bench_malloc(maxdim * sizeof(bench_real));
     HARM_F77(a, m, inv, s, &ifset, &iferr);
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, -1);
}

void doit(int iter, struct problem *p)
{
     int i, ifset = p->sign < 0 ? -2 : 2, iferr = 0;
     bench_complex *a = (bench_complex *) p->in;

     for (i = 0; i < iter; ++i) {
	  HARM_F77(a, m, inv, s, &ifset, &iferr);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(s);
     bench_free(inv);
}

