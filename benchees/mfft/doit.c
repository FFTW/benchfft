/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "mfft")
BENCH_DOC("author", "A. Nobile and V. Roberto")
BENCH_DOC("year", "1987")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("url", "http://risc2.numis.nwu.edu/ftp/pub/transforms/mfft.tar.gz")
BENCH_DOC("url-was-valid-on", "Sat Aug 31 16:35:13 EDT 2002")
BENCH_DOC("bibitem", 
	  "A. Nobile and V. Roberto, MFFT: A package for two- and three-dimensional vectorized discrete Fourier transforms, Computer Physics Communications 42, p. 233 (1986).")
END_BENCH_DOC

#define C2FFT F77_FUNC(c2fft,C2FFT)
extern void C2FFT(bench_complex *x, unsigned int *ldn, 
		  unsigned int *n1, unsigned int *n2,
		  int *w1, int *w2,
		  int *isig, int *iord, int *iwork, int *ierr);

#define C3FFT F77_FUNC(c3fft,C3FFT)
extern void C3FFT(bench_complex *x, unsigned int *ldn, 
		  unsigned int *n1, unsigned int *n2, unsigned int *n3,
		  int *w1, int *w2, int *w3,
		  int *iopt, int *isig, int *iord, int *iwork, int *ierr);

#define C4FFT F77_FUNC(c4fft,C4FFT)
extern void C4FFT(bench_complex *x, 
		  unsigned int *ldn1, unsigned int *ldn2,
		  unsigned int *n1, unsigned int *n2, unsigned int *n3,
		  unsigned int *n4,
		  int *w1, int *w2, int *w3, int *w4,
		  int *isig, int *iord, int *iwork, int *ierr);

#define R2FFT F77_FUNC(r2fft,R2FFT)
extern void R2FFT(bench_complex *x, unsigned int *ldn, 
		  unsigned int *n1, unsigned int *n2,
		  int *w1, int *w2,
		  int *isig, int *iord, int *iwork, int *ierr);

#define R3FFT F77_FUNC(r3fft,R3FFT)
extern void R3FFT(bench_complex *x, unsigned int *ldn, 
		  unsigned int *n1, unsigned int *n2, unsigned int *n3,
		  int *w1, int *w2, int *w3,
		  int *iopt, int *isig, int *iord, int *iwork, int *ierr);

#define R4FFT F77_FUNC(r4fft,R4FFT)
extern void R4FFT(bench_complex *x, 
		  unsigned int *ldn1, unsigned int *ldn2,
		  unsigned int *n1, unsigned int *n2, unsigned int *n3,
		  unsigned int *n4,
		  int *w1, int *w2, int *w3, int *w4,
		  int *isig, int *iord, int *iwork, int *ierr);

int n_ok(unsigned int rank, unsigned int *n)
{
     unsigned int i;
     for (i = 0; i < rank; ++i)
	  if (!check_prime_factors(n[i], 5))
	       return 0;
     return 1;
}

int can_do(struct problem *p)
{
     return (SINGLE_PRECISION &&
	     sizeof(float) == sizeof(int) && /* thanks to ugly hack */
	     p->rank >= 2 && p->rank <= 4 &&
	     (p->kind == PROBLEM_COMPLEX || p->n[p->rank - 1] % 2 == 0) &&
	     problem_in_place(p) && n_ok(p->rank, p->n));
}

#define IOPT 0

static unsigned int ldn1 = 1, n1 = 1, n2 = 1, n3 = 1, n4 = 1;
static int *w1, *w2, *w3, *w4, *iwork;
static int iord = 1, ierr, iopt = IOPT;

#define MAX2(a,b) ((a) > (b) ? (a) : (b))

void setup(struct problem *p)
{
     int isig = 0; /* indicates initialization call */

     BENCH_ASSERT(can_do(p));

     switch (p->rank) {
	 case 2: n2 = p->n[0]; n1 = p->n[1]; break;
	 case 3: n3 = p->n[0]; n2 = p->n[1]; n1 = p->n[2]; break;
	 case 4: n4 = p->n[0]; n3 = p->n[1]; n2 = p->n[2]; n1 = p->n[3]; break;
     }
     if (p->kind == PROBLEM_COMPLEX)
	  ldn1 = n1;
     else
	  ldn1 = n1 + 2;

     w1 = (int *) bench_malloc(sizeof(int) *
			       ((p->kind == PROBLEM_REAL ? 6 : 4) * n1 + 14));
     if (iopt == 1 && p->rank == 3)
	  w2 = (int *) bench_malloc(sizeof(int) * (4*n2*(n1+1) + 14));
     else
	  w2 = (int *) bench_malloc(sizeof(int) * (4*n2 + 14));
     w3 = (int *) bench_malloc(sizeof(int) * (4*n3 + 14));
     w4 = (int *) bench_malloc(sizeof(int) * (4*n4 + 14));
     iwork = (int *) bench_malloc(sizeof(int) * MAX2(n1,MAX2(n2,MAX2(n3,n4))));
     
     if (p->kind == PROBLEM_COMPLEX) {
	  switch (p->rank) {
	      case 2:
		   C2FFT((bench_complex*) p->in, &ldn1, &n1, &n2,
			 w1, w2, &isig, &iord, iwork, &ierr);
		   break;
	      case 3:
		   C3FFT((bench_complex*) p->in, &ldn1, &n1, &n2, &n3,
			 w1, w2, w3, &iopt, &isig, &iord, iwork, &ierr);
		   break;
	      case 4:
		   C4FFT((bench_complex*) p->in,
			 &ldn1, &n2, &n1, &n2, &n3, &n4,
			 w1, w2, w3, w4, &isig, &iord, iwork, &ierr);
		   break;
	  }
     }
     else /* PROBLEM_REAL */ {
	  switch (p->rank) {
	      case 2:
		   R2FFT((bench_complex*) p->in, &ldn1, &n1, &n2,
			 w1, w2, &isig, &iord, iwork, &ierr);
		   break;
	      case 3:
		   R3FFT((bench_complex*) p->in, &ldn1, &n1, &n2, &n3,
			 w1, w2, w3, &iopt, &isig, &iord, iwork, &ierr);
		   break;
	      case 4:
		   R4FFT((bench_complex*) p->in,
			 &ldn1, &n2, &n1, &n2, &n3, &n4,
			 w1, w2, w3, w4, &isig, &iord, iwork, &ierr);
		   break;
	  }
     }
     BENCH_ASSERT(ierr == 0);
}

void doit(int iter, struct problem *p)
{
     bench_complex *in = p->in;
     int isig = p->sign;
     int i;

     if (p->kind == PROBLEM_COMPLEX) {
	  switch (p->rank) {
	      case 2:
		   for (i = 0; i < iter; ++i)
			C2FFT(in, &ldn1, &n1, &n2,
			      w1, w2, &isig, &iord, iwork, &ierr);
		   break;
	      case 3:
		   for (i = 0; i < iter; ++i)
			C3FFT(in, &ldn1, &n1, &n2, &n3,
			      w1, w2, w3, &iopt, &isig, &iord, iwork, &ierr);
		   break;
	      case 4:
		   for (i = 0; i < iter; ++i)
			C4FFT(in, &ldn1, &n2, &n1, &n2, &n3, &n4,
			      w1, w2, w3, w4, &isig, &iord, iwork, &ierr);
		   break;
	  }
     }
     else /* PROBLEM_REAL */ {
	  switch (p->rank) {
	      case 2:
		   for (i = 0; i < iter; ++i)
			R2FFT(in, &ldn1, &n1, &n2,
			      w1, w2, &isig, &iord, iwork, &ierr);
		   break;
	      case 3:
		   for (i = 0; i < iter; ++i)
			R3FFT(in, &ldn1, &n1, &n2, &n3,
			      w1, w2, w3, &iopt, &isig, &iord, iwork, &ierr);
		   break;
	      case 4:
		   for (i = 0; i < iter; ++i)
			R4FFT(in, &ldn1, &n2, &n1, &n2, &n3, &n4,
			      w1, w2, w3, w4, &isig, &iord, iwork, &ierr);
		   break;
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     free(iwork);
     free(w4);
     free(w3);
     free(w2);
     free(w1);
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_unpacked(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_unpacked(p, in, -1.0);
}

void copy_r2c(struct problem *p, bench_complex *out)
{
     copy_r2c_unpacked(p, out);
}

void copy_c2r(struct problem *p, bench_complex *in)
{
     copy_c2r_unpacked(p, in);
}
