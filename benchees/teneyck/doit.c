/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "teneyck")
BENCH_DOC("author", "Lynn F. Ten Eyck")
BENCH_DOC("year", "1973")
BENCH_DOC("email", "lteneyck@sdsc.edu")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("url", "http://risc2.numis.nwu.edu/ftp/pub/transforms/fftlib.f")
BENCH_DOC("url-was-valid-on", "Sat Aug 31 16:37:44 EDT 2002")
BENCH_DOC("bibitem",
	  "Modified by L. F. Ten Eyck from a one-dimensional version written by G. T. Sande, 1969.")
BENCH_DOC("bibitem", 
	  "L. F. Ten Eyck,  Crystallographic Fast Fourier Transforms, Acta Cryst. A29, 183-191 (1973).")
END_BENCH_DOC

int n_ok(unsigned int rank, unsigned int *n)
{
     unsigned int i;
     for (i = 0; i < rank; ++i) {
	  int num_facts = check_prime_factors(n[i], 19);
	  if (!num_facts || num_facts > 14)
               return 0;
     }
     return 1;
}

int can_do(struct problem *p)
{
     return (SINGLE_PRECISION &&
	     p->rank >= 1 && p->rank <= 3 &&
	     p->kind == PROBLEM_COMPLEX &&
	     n_ok(p->rank, p->n) &&
	     problem_in_place(p)
	  );
}

#define CMPLFT F77_FUNC(cmplft,CMPLFT)

extern void CMPLFT(bench_real *x, bench_real *y, int *n, int d[5]);

int d0[5] = {0,0,0,0,0}, d1[5] = {0,0,0,0,0}, d2[5] = {0,0,0,0,0};

/* initialize CMPLFT's screwy dims parameter */
void set_dims(int d[5], int nbefore, int n, int nafter, int stride)
{
     int i;
     d[0] = nbefore * n * nafter;
     d[1] = nafter;
     d[2] = n * nafter;
     d[3] = nafter;
     d[4] = 1;
     for (i = 0; i < 5; ++i)
	  d[i] *= stride;
}

void copy_c2c_from(struct problem *p, bench_complex *in)
{
     if (p->split) {
	  bench_real *x = (bench_real *) p->in;
	  copy_c2ri(in, x, x + p->size, p->size);
     }
     else
	  cacopy(in, p->in, p->size);
}

void copy_c2c_to(struct problem *p, bench_complex *out)
{
     if (p->split) {
	  bench_real *x = (bench_real *) p->out;
	  copy_ri2c(x, x + p->size, out, p->size);
     }
     else
	  cacopy(p->out, out, p->size);
}

void setup(struct problem *p)
{
     int stride;
     if (p->split)
	  stride = 1;
     else
	  stride = 2;
     switch (p->rank) {
	 case 1:
	      set_dims(d0, 1, p->n[0], 1, stride);
	      break;
	 case 2:
	      set_dims(d0, 1, p->n[0], p->n[1], stride);
	      set_dims(d1, p->n[0], p->n[1], 1, stride);
	      break;
	 case 3:
	      set_dims(d0, 1, p->n[0], p->n[1] * p->n[2], stride);
	      set_dims(d1, p->n[0], p->n[1], p->n[2], stride);
	      set_dims(d2, p->n[0] * p->n[1], p->n[2], 1, stride);
	      break;
	 default:
	      BENCH_ASSERT(0);
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     bench_real *x, *y;

     if (p->sign < 0) {
	  x = (bench_real *) p->in;
	  y = x + (p->split ? p->size : 1);
     }
     else {  /* inverse transform is given by swapping real/imag parts */
	  y = (bench_real *) p->in;
	  x = y + (p->split ? p->size : 1);
     }

     switch (p->rank) {
	 case 1: {
	      int n0 = p->n[0];
	      for (i = 0; i < iter; ++i) {
		   CMPLFT(x, y, &n0, d0);
	      }
	      break;
	 }
	 case 2: {
	      int n0 = p->n[0], n1 = p->n[1];
	      for (i = 0; i < iter; ++i) {
		   CMPLFT(x, y, &n0, d0);
		   CMPLFT(x, y, &n1, d1);
	      }
	      break;
	 }
	 case 3: {
	      int n0 = p->n[0], n1 = p->n[1], n2 = p->n[2];
	      for (i = 0; i < iter; ++i) {
		   CMPLFT(x, y, &n0, d0);
		   CMPLFT(x, y, &n1, d1);
		   CMPLFT(x, y, &n2, d2);
	      }
	      break;
	 }
	 default:
	      BENCH_ASSERT(0);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}

