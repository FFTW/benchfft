/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC

BENCH_DOC("package", "Numerical Recipes")
#ifdef FORTRAN
  BENCH_DOC("name", "nr-f")
  BENCH_DOC("language", "Fortran 77")
  BENCH_DOC("year", "1992")
  BENCH_DOC("bibitem", 
	    "W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery, Numerical Recipes in Fortran (Cambridge Univ. Press, 1992).")
#else
  BENCH_DOC("name", "nr-c")
  BENCH_DOC("language", "C")
  BENCH_DOC("year", "1993")
  BENCH_DOC("bibitem", 
	    "W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery, Numerical Recipes in C: The Art of Scientific Computing (Cambridge Univ. Press, 1993).")
#endif

BENCH_DOC("author", "W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery")
BENCH_DOC("url", "http://www.nr.com/")
BENCH_DOC("url-was-valid-on", "Fri Aug 17 15:35:44 EDT 2001")

BENCH_DOC("copyright", "Available under a restrictive license, but may not be redistributed in source form.")
END_BENCH_DOC


int can_do(struct problem *p)
{
     return (SINGLE_PRECISION &&
	     p->rank > 0 &&
	     (p->kind == PROBLEM_COMPLEX || 
	      p->n[p->rank - 1] % 2 == 0) &&
	     (p->kind == PROBLEM_COMPLEX || p->rank <= 3) &&
	     problem_power_of_two(p, 1)
	  );
}

#ifdef FORTRAN

typedef int size_int;

#define FOUR1 F77_FUNC(four1,FOUR1)
extern void FOUR1(float *data, size_int *nn, int *isign);
#define fft(data,n,is) FOUR1(data,&n,&is)

#define FOURN F77_FUNC(fourn,FOURN)
extern void FOURN(float *data, size_int *nn, int *ndim, int *isign);
#define fftn(data,nn,rank,is) FOURN(data,nn,&rank,&is)

#define REALFT F77_FUNC(realft,REALFT)
extern void REALFT(float *data, size_int *nn, int *isign);
#define rfft(data,n,is) REALFT(data,&n,&is)

#define RLFT3 F77_FUNC(rlft3,RLFT3)
extern void RLFT3(float *data, float *speq, 
		  size_int *nn1, size_int *nn2, size_int *nn3,
		  int *isign);
#define rfft3(data,speq,n1,n2,n3,is) RLFT3(data,speq,&n3,&n2,&n1,&is)

#else /* C */

typedef unsigned long size_int;

extern void four1(float *data, size_int nn, int isign);
#define fft(data,n,is) four1(data-1,n,is)

extern void fourn(float *data, size_int *nn, int ndim, int isign);
#define fftn(data,nn,rank,is) fourn(data-1,nn-1,rank,is)

extern void realft(float *data, size_int nn, int isign);
#define rfft(data,n,is) realft(data-1,n,is)

extern void rlft3(float ***data, float **speq, 
		  size_int nn1, size_int nn2, size_int nn3,
		  int isign);
#define rfft3(data,speq,n1,n2,n3,is) rlft3(data##_p,speq##_p,n1,n2,n3,is)

#endif

size_int n, n1 = 1, n2 = 1, n3 = 1, nn[30];
float *speq = 0;
#ifndef FORTRAN
float ***data_p, **speq_p;
#endif

void setup(struct problem *p)
{
     unsigned int i;
     BENCH_ASSERT(can_do(p));
     
     n = p->size;

#ifdef FORTRAN
     for (i = 0; i < p->rank; ++i)
	  nn[p->rank - 1 - i] = p->n[i];
#else
     for (i = 0; i < p->rank; ++i)
	  nn[i] = p->n[i];
#endif
     if (p->rank == 3) {
	  n1 = p->n[0]; n2 = p->n[1]; n3 = p->n[2];
     }
     else if (p->rank == 2) {
	  n1 = 1; n2 = p->n[0]; n3 = p->n[1];
     }

     if (p->kind == PROBLEM_REAL && p->rank > 1) {
	  float *data = (float *) p->in;
	  speq = (float *) bench_malloc(sizeof(float) * n1 * (n2 * 2));
#ifndef FORTRAN
	  speq_p = (float **) bench_malloc(sizeof(float*) * (n1 + 1));
	  data_p = (float ***) bench_malloc(sizeof(float**) * (n1 + 1));
	  for (i = 1; i <= n1; ++i) {
	       unsigned int j;
	       data_p[i] = (float **) bench_malloc(sizeof(float*) * (n2+1));
	       for (j = 1; j <= n2; ++j)
		    data_p[i][j] = data + ((i-1)*n2 + (j-1))*n3 - 1;
	       speq_p[i] = speq + (i-1)*n2*2 - 1;
	  }
#endif
     }
}

void done(struct problem *p)
{
     unsigned int i;
     if (p->kind == PROBLEM_REAL && p->rank > 1) {
	  bench_free(speq);
#ifndef FORTRAN
	  for (i = 1; i <= n1; ++i)
	       bench_free(data_p[i]);
	  bench_free(data_p);
	  bench_free(speq_p);
#endif
     }
}

void doit(int iter, struct problem *p)
{
     float *data = (float *) p->in;
     int i;

     if (p->kind == PROBLEM_COMPLEX) {
	  int isign = p->sign;
	  if (p->rank == 1)
	       for (i = 0; i < iter; ++i) {
		    fft(data, n, isign);
	       }
	  else {
	       int rank = p->rank;
	       for (i = 0; i < iter; ++i) {
		    fftn(data, nn, rank, isign);
	       }
	  }
     }
     else {
	  int isign = -p->sign;
	  if (p->rank == 1)
	       for (i = 0; i < iter; ++i) {
		    rfft(data, n, isign);
	       }
	  else {
	       for (i = 0; i < iter; ++i) {
		    rfft3(data, speq, n1, n2, n3, isign);
	       }
	  }
     }
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     if (p->rank == 1) {
	  copy_h2c_1d_packed(p, out, 1.0);
     }
     else {
	  int i, j, k;
	  int n0, n1, n2;
	  bench_complex *cdata = (bench_complex *) p->out;
	  bench_complex *cspeq = (bench_complex *) speq;

	  if (p->rank == 3) {
	       n0 = p->n[0]; n1 = p->n[1]; n2 = p->n[2];
	  }
	  else {
	       n0 = 1; n1 = p->n[0]; n2 = p->n[1];
	  }

	  for (i = n0 - 1; i >= 0; --i)
	       for (j = n1 - 1; j >= 0; --j) {
		    for (k = n2/2 - 1; k >= 0; --k) {
			 cdata[(i*n1 + j)*(n2/2+1) + k] =
			      cdata[(i*n1 + j)*(n2/2) + k];
		    }
		    cdata[(i*n1 + j)*(n2/2+1) + n2/2] =
			 cspeq[i*n1 + j];
	       }

	  copy_h2c_unpacked(p, out, 1.0);
     }
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     if (p->rank == 1) {
	  copy_c2h_1d_packed(p, in, 1.0);
     }
     else {
	  int i, j, k;
	  int n0, n1, n2;
	  bench_complex *cdata = (bench_complex *) p->in;
	  bench_complex *cspeq = (bench_complex *) speq;

	  if (p->rank == 3) {
	       n0 = p->n[0]; n1 = p->n[1]; n2 = p->n[2];
	  }
	  else {
	       n0 = 1; n1 = p->n[0]; n2 = p->n[1];
	  }

	  for (i = 0; i < n0; ++i)
	       for (j = 0; j < n1; ++j) {
		    for (k = 0; k < n2/2; ++k) {
			 cdata[(i*n1 + j)*(n2/2) + k] = in[(i*n1 + j)*n2 + k];
			 c_im(cdata[(i*n1 + j)*(n2/2) + k]) *= -1;
		    }
		    cspeq[i*n1 + j] = in[(i*n1 + j)*n2 + k];
		    c_im(cspeq[i*n1 + j]) *= -1;
	       }
     }
}

/* Backwards real transform is scaled by 0.5 */
void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     if (p->kind == PROBLEM_REAL && p->sign > 0) {
          bench_complex x;
          c_re(x) = 2.0;
          c_im(x) = 0;

          cascale(out, p->size, x);
     }
}
