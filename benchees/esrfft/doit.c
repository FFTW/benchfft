/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "esrfft")
BENCH_DOC("author", "Daisuke Takahashi")
BENCH_DOC("year", "2002")
BENCH_DOC("email", "daisuke@is.tsukuba.ac.jp")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("url", "http://www.ffte.jp/")
BENCH_DOC("url-was-valid-on", "Fri Jan  3 16:11:54 EST 2003")
BENCH_DOC("notes", "Complex data are stored in separate real/imag arrays.")
BENCH_DOC("bibitem", 
	  "D. Takahashi, An extended split-radix FFT algorithm, IEEE Signal Processing Lett. 8, pp. 145-147, May 2001.")
BENCH_DOC("notes", "The backward transform is scaled, and also has an additional pre-processing pass to flip the sign of the input imaginary part.  Since the code accepts the real/imag parts as separate arrays, we instead swap the arguments to avoid this overhead.")
BENCH_DOC("copyright",
	  "Copyright(C) 2000-2002 Daisuke Takahashi (e-mail: daisuke@is.tsukuba.ac.jp or ffte@ffte.jp)\n\n"
	  "You may use, copy, modify this code for any purpose (include commercial use) and without fee. You may distribute this ORIGINAL package.")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (DOUBLE_PRECISION &&
	     p->rank == 1 &&
	     p->kind == PROBLEM_COMPLEX &&
	     problem_power_of_two(p, 1)
	  );
}

#define ESRFFT F77_FUNC(esrfft,ESRFFT)

extern void ESRFFT(bench_real *x, bench_real *y, 
		   int *ibeta, int *n, int *m, int *isign);

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

int m, *ibeta;

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
     m = log_2(p->n[0]);
     ibeta = (int *) bench_malloc((p->n[0] / 2) * sizeof(int));
}

void doit(int iter, struct problem *p)
{
     int plusone = 1;
     int i;
     int n = p->n[0];
     bench_real *x, *y;

     if (p->sign < 0) {
	  x = (bench_real *) p->in;
	  y = x + n;
     }
     else {  /* inverse transform is given by swapping real/imag parts */
	  y = (bench_real *) p->in;
	  x = y + n;
     }
     
     for (i = 0; i < iter; ++i) {
	  ESRFFT(x, y, ibeta, &n, &m, &plusone);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(ibeta);
}
