/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "as117")
BENCH_DOC("package", "StatLib")
BENCH_DOC("author", "Donald M. Monro and John L. Branch")
BENCH_DOC("year", "1977")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("url", "http://lib.stat.cmu.edu/apstat/117")
BENCH_DOC("url-was-valid-on", "Sun Jan 13 14:14:59 EST 2002")
BENCH_DOC("notes", "Downloaded from the StatLib repository at Carnegie Mellon University, in the Applied Statistics algorithm collection.")
BENCH_DOC("notes", "Complex data are stored in separate real/imag arrays.")
BENCH_DOC("notes", "Uses the chirp-z algorithm to compute FFTs of arbitrary size, via the as83 power-of-two FFT.")
BENCH_DOC("notes", "The forward transform is scaled.")
BENCH_DOC("bibitem", 
	  "Donald M. Monro and John L. Branch, Algorithm AS 117: The chirp discrete Fourier transform of general length, Appl. Statist. 26 (3), pp. 351-361 (1977).")
BENCH_DOC("copyright", "The Royal Statistical Society holds the copyright to these routines, but has given its permission for their distribution provided that no fee is charged.")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (SINGLE_PRECISION &&
	     p->rank == 1 &&
	     ((p->n[0] >= 3 && p->n[0] <= (1<<19) - 1 &&
	      p->kind == PROBLEM_COMPLEX) ||
	      (p->n[0] >= 6 && p->n[0] <= (1<<20) - 2 &&
	       p->n[0] % 2 == 0 &&
	       p->kind == PROBLEM_REAL)) &&
	     problem_in_place(p)
	  );
}

#define CHRFT F77_FUNC(chrft,CHRFT)
extern void CHRFT(bench_real *re, bench_real *im,
		  bench_real *wre, bench_real *wim,
		  unsigned int *n, unsigned int *nwork, int *isn, int *ierr);

#define CHFOR F77_FUNC(chfor,CHFOR)
extern void CHFOR(bench_real *x,
		  bench_real *wre, bench_real *wim,
		  unsigned int *n, unsigned int *nwork2, unsigned int *nwork,
		  int *ierr);

#define CHREV F77_FUNC(chrev,CHREV)
extern void CHREV(bench_real *x,
		  bench_real *wre, bench_real *wim,
		  unsigned int *n, unsigned int *nwork2, unsigned int *nwork,
		  int *ierr);

#define SETWT F77_FUNC(setwt,SETWT)
extern void SETWT(bench_real *wre, bench_real *wim,
		  unsigned int *n, unsigned int *nwork, int *ierr);

unsigned int nwork = 0;
bench_real *xre = 0, *xim = 0;
bench_real *wre = 0, *wim = 0;

void setup(struct problem *p)
{
     unsigned int n;
     BENCH_ASSERT(can_do(p));

     n = p->n[0];
     if (p->kind == PROBLEM_REAL)
	  n /= 2;
     nwork = 1;
     while (nwork < 2 * n)
	  nwork *= 2;
     wre = bench_malloc(nwork * sizeof(bench_real));
     wim = bench_malloc(nwork * sizeof(bench_real));
     xre = bench_malloc(nwork * sizeof(bench_real) * 2);
     xim = xre + nwork;
     {
	  int ierr;
	  SETWT(wre, wim, &n, &nwork, &ierr);
	  BENCH_ASSERT(ierr == 0);
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     unsigned int n = p->n[0];
     int ierr;

     if (p->kind == PROBLEM_COMPLEX) {
	  int isn = (p->sign > 0) ? -1 : 1;
	  for (i = 0; i < iter; ++i) {
	       CHRFT(xre, xim, wre, wim, &n, &nwork, &isn, &ierr);
	  }
     }
     else if (p->sign < 0) { /* forward real transform */
	  unsigned int nwork2 = nwork * 2;
          for (i = 0; i < iter; ++i) {
               CHFOR(xre, wre, wim, &n, &nwork2, &nwork, &ierr);
          }
     }
     else { /* backward real transform */
	  unsigned int nwork2 = nwork * 2;
          for (i = 0; i < iter; ++i) {
               CHREV(xre, wre, wim, &n, &nwork2, &nwork, &ierr);
          }
     }
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, -1);
}

void copy_c2c_from(struct problem *p, bench_complex *in)
{
     copy_c2ri(in, xre, xim, p->size);
}

void copy_c2c_to(struct problem *p, bench_complex *out)
{
     copy_ri2c(xre, xim, out, p->size);
}

void copy_c2r(struct problem *p, bench_complex *in)
{
     unsigned int i;
     unsigned int n = p->size;

     for (i = 0; i < n; ++i)
          xre[i] = c_re(in[i]);
}

void copy_r2c(struct problem *p, bench_complex *out)
{
     unsigned int i;
     unsigned int n = p->size;

     for (i = 0; i < n; ++i) {
          c_re(out[i]) = xre[i];
          c_im(out[i]) = 0.0;
     }
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     int i, n = p->size;
     bench_real *x = (bench_real *) p->out;

     for (i = 0; i < n; ++i)
	  x[i] = xre[i];
     for (i = n/2 + 1; i < n - (i - n/2); ++i) {
          bench_real y = x[i];
          x[i] = x[n - (i - n/2)];
          x[n - (i - n/2)] = y;
     }

     copy_h2c_1d_halfcomplex(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     int i, n = p->size;
     bench_real *x = (bench_real *) p->in;

     copy_c2h_1d_halfcomplex(p, in, -1.0);

     for (i = n/2 + 1; i < n - (i - n/2); ++i) {
          bench_real y = x[i];
          x[i] = x[n - (i - n/2)];
          x[n - (i - n/2)] = y;
     }
     for (i = 0; i < n; ++i)
	  xre[i] = x[i];
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(xre);
     bench_free(wim);
     bench_free(wre);
}
