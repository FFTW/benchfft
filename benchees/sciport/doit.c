/* this program is in the public domain */

#include "bench-user.h"

BEGIN_BENCH_DOC
BENCH_DOC("name", "sciport")
BENCH_DOC("author", "Scott H. Lamson")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("url", "http://www.netlib.org/scilib/")
BENCH_DOC("url-was-valid-on", "Fri Aug 30 21:09:06 EDT 2002")
BENCH_DOC("bibitem", 
	  "SCIPORT is a portable re-implementation of Cray's SCILIB; "
	  "developed at General Electric, probably by Scott H. Lamson.")
BENCH_DOC("notes", "Stockham auto-sort FFT, radix-2.")
BENCH_DOC("notes", "Slightly modified to compile properly with g77.")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (p->rank == 1 &&
	     p->n[0] >= 4 &&
	     (p->kind == PROBLEM_COMPLEX || p->n[0] >= 8) &&
	     power_of_two(p->n[0])
	  );
}

#define CFFT2 F77_FUNC(cfft2,CFFT2)
#define RCFFT2 F77_FUNC(rcfft2,RCFFT2)
#define CRFFT2 F77_FUNC(crfft2,CRFFT2)

typedef unsigned int size_int;
extern void CFFT2(int *init, int *isign, size_int *n, bench_complex *cx,
		  bench_complex *cwork, bench_complex *cy);
extern void RCFFT2(int *init, int *isign, size_int *n, bench_complex *cx,
		   bench_complex *cwork, bench_complex *cy);
extern void CRFFT2(int *init, int *isign, size_int *n, bench_complex *cx,
		   bench_complex *cwork, bench_complex *cy);

bench_complex *work = 0;

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_unpacked(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_unpacked(p, in, -1.0);
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     if (p->kind == PROBLEM_REAL && p->sign < 0) {
	  bench_complex x;
	  c_re(x) = 0.5;
	  c_im(x) = 0;
	  cascale(out, p->size, x);
     }
}

void setup(struct problem *p)
{
     int init = 1, isign = -1;
     size_int n;
     bench_complex *a = (bench_complex *) p->in;
     bench_complex *b = (bench_complex *) p->out;
     BENCH_ASSERT(can_do(p));

     n = p->n[0];
     if (p->kind == PROBLEM_COMPLEX) {
	  work = (bench_complex*) bench_malloc((2*n + n/2) 
					       * sizeof(bench_complex));
	  CFFT2(&init, &isign, &n, a, work, b);
     }
     else {
	  work = (bench_complex*) bench_malloc((n + 2 + n/2) 
					       * sizeof(bench_complex));
	  if (p->sign < 0)
	       RCFFT2(&init, &isign, &n, a, work, b);
	  else
	       CRFFT2(&init, &isign, &n, a, work, b);
     }
}

void doit(int iter, struct problem *p)
{
     int i, init = 0, isign = p->sign;
     size_int n = p->n[0];
     bench_complex *a = (bench_complex *) p->in;
     bench_complex *b = (bench_complex *) p->out;

     if (p->kind == PROBLEM_COMPLEX)
	  for (i = 0; i < iter; ++i)
	       CFFT2(&init, &isign, &n, a, work, b);
     else if (isign < 0)
	  for (i = 0; i < iter; ++i)
	       RCFFT2(&init, &isign, &n, a, work, b);
     else
	  for (i = 0; i < iter; ++i)
	       CRFFT2(&init, &isign, &n, a, work, b);
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(work);
}

