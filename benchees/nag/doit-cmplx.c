/* this program is in the public domain */

#include "bench-user.h"

#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "nag-cmplx")
BENCH_DOC("package", "NAG Fortran Library")
BENCH_DOC("author", "Numerical Algorithms Group")
BENCH_DOC("url", "http://www.nag.com/")
BENCH_DOC("url-was-valid-on", "Wed Aug 15 13:43:43 EDT 2001")
BENCH_DOC("notes", "Uses c06pcf and friends: complex data type (interleaved real/imaginary parts), extra workspace supplied.")
BENCH_DOC("notes", "both forward and backward transforms are scaled")
BENCH_DOC("copyright", "Proprietary, commercial software.")
END_BENCH_DOC

static int ifail = 0;

/* Don't you love the beautiful 6-character Fortran-77 identifiers? */

#define C06PCF F77_FUNC(c06pcf,C06PCF)
extern void C06PCF(char *direct, bench_complex *x, unsigned int *n,
                   double *work, int *ifail);
#define fft(d,n,x,work) C06PCF(&d, x, &n, work, &ifail)

#define C06PAF F77_FUNC(c06paf,C06PAF)
extern void C06PAF(char *direct, double *x, unsigned int *n, double *work, int *ifail);
#define fftrc(d,n,x,work) C06PAF(&d, x, &n, work, &ifail)

#define C06PUF F77_FUNC(c06puf,C06PUF)
extern void C06PUF(char *direct, unsigned int *m, unsigned int *n, 
		   bench_complex *x, double *work, int *ifail);
#define fft2(d,n0,n1,x,w) C06PUF(&d, &n1,&n0, x, w, &ifail)

#define C06PXF F77_FUNC(c06pxf,C06PXF)
extern void C06PXF(char *direct,
		   unsigned int *n1, unsigned int *n2, unsigned int *n3,
		   bench_complex *x, double *work, int *ifail);
#define fft3(d,n0,n1,n2,x,w) C06PXF(&d, &n2,&n1,&n0, x, w, &ifail)

int n_ok(unsigned int rank, unsigned int *n)
{
     int i;
     for (i = 0; i < rank; ++i)
	  if (n[i] <= 1)
	       return 0;
     return 1;
}

int can_do(struct problem *p)
{
     return (DOUBLE_PRECISION &&
	     p->rank >= 1 &&
	     ((p->rank <= 3 && p->kind == PROBLEM_COMPLEX) || (p->rank == 1))
	     && problem_in_place(p) && n_ok(p->rank, p->n));
}

double *work = 0;

void setup(struct problem *p)
{
     int nwork = 0;
     BENCH_ASSERT(can_do(p));

     if (p->kind == PROBLEM_REAL) {
	  nwork = 2*p->size + 15;
     }
     else {
	  nwork = 15*p->rank + p->size + p->n[0];
	  if (p->rank >= 2)
	       nwork += p->n[1];
	  if (p->rank >= 3)
	       nwork += p->n[2];
	  nwork *= 2; /* 1 complex == 2 double */
     }
     work = (double*) bench_malloc(sizeof(double) * nwork);
}

void doit(int iter, struct problem *p)
{
     int i;
     char d = p->sign < 0 ? 'F' : 'B';

     if (p->kind == PROBLEM_COMPLEX) {
	  bench_complex *x = (bench_complex *) p->in;
	  
	  switch (p->rank) {
	      case 1:
	      {
		   unsigned int n0 = p->n[0];
		   for (i = 0; i < iter; ++i) {
			fft(d, n0, x, work);
		   }
		   break;
	      }
	      case 2:
	      {
		   unsigned int n0 = p->n[0], n1 = p->n[1];
		   for (i = 0; i < iter; ++i) {
			fft2(d, n0, n1, x, work);
		   }
		   break;
	      }
	      case 3:
	      {
		   unsigned int n0 = p->n[0], n1 = p->n[1], n2 = p->n[2];
		   for (i = 0; i < iter; ++i) {
			fft3(d, n0, n1, n2, x, work);
		   }
		   break;
	      }
	  }
     }
     else /* PROBLEM_REAL */ {
	  double *x = (double *) p->in;
	  unsigned int n0 = p->n[0];
	  
	  for (i = 0; i < iter; ++i) {
	       fftrc(d, n0, x, work);
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(work);
}

/* Undo normalization: */
void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     bench_complex x;
     c_re(x) = sqrt(p->size);
     c_im(x) = 0;
     
     cascale(out, p->size, x);
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_unpacked(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_unpacked(p, in, -1.0);
}
