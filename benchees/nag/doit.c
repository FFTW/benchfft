/* this program is in the public domain */

#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", NAME)
BENCH_DOC("package", "NAG Fortran Library")
BENCH_DOC("author", "Numerical Algorithms Group")
BENCH_DOC("url", "http://www.nag.com/")
BENCH_DOC("url-was-valid-on", "Wed Aug 15 13:43:43 EDT 2001")
BENCH_DOC("notes", NOTE)
BENCH_DOC("notes", "Complex data are stored in separate real/imag arrays.")
BENCH_DOC("notes", "both forward and backward transforms are scaled")
BENCH_DOC("notes", "forwards and backwards real transforms have the same sign")
BENCH_DOC("copyright", "Proprietary, commercial software.")
END_BENCH_DOC

static unsigned int one = 1;
static int ifail = 0;
static char INIT = 'I', SUBS = 'S';

/* Don't you love the beautiful 6-character Fortran-77 identifiers? */

#define C06FUF F77_FUNC(c06fuf,C06FUF)
extern void C06FUF(unsigned int *m, unsigned int *n, double *x, double *y, char *init,
		   double *trigm, double *trign, double *work, int *ifail);
#define fft2(n0,n1,x,y,t0,t1,w) C06FUF(&n1,&n0, x, y, &SUBS, t1,t0, w, &ifail)
#define init2(n0,n1,x,y,t0,t1,w) C06FUF(&n1,&n0, x, y, &INIT, t1,t0, w, &ifail)

#define C06FXF F77_FUNC(c06fxf,C06FXF)
extern void C06FXF(unsigned int *n1, unsigned int *n2, unsigned int *n3, double *x, double *y, char *init,
		   double *trig1, double *trig2, double *trig3,
		   double *work, int *ifail);
#define fft3(n0,n1,n2,x,y,t0,t1,t2,w) C06FXF(&n2,&n1,&n0, x, y, &SUBS, t2,t1,t0, w, &ifail)
#define init3(n0,n1,n2,x,y,t0,t1,t2,w) C06FXF(&n2,&n1,&n0, x, y, &INIT, t2,t1,t0, w, &ifail)

#define C06GQF F77_FUNC(c06gqf,C06GQF)
extern void C06GQF(unsigned int *m, unsigned int *n, double *x, int *ifail);
#define conj_hc(n,x) C06GQF(&one, &n, x, &ifail)

int n_ok(unsigned int rank, unsigned int *n)
{
     int i;
     for (i = 0; i < rank; ++i)
	  if (n[i] <= 1 || !check_prime_factors(n[i], 19))
	       return 0;
     return 1;
}

int can_do(struct problem *p)
{
     return (DOUBLE_PRECISION &&
	     p->rank >= 1 &&
	     ((p->rank <= 3 && p->kind == PROBLEM_COMPLEX) || (p->rank == 1))
#ifdef ONLY_1D
	     && p->rank == 1
#endif
	     && problem_in_place(p) && n_ok(p->rank, p->n));
}

double *trig0 = 0, *trig1 = 0, *trig2 = 0;
double *work = 0;

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));

     if (p->kind == PROBLEM_REAL) {
	  work = (double*) bench_malloc(sizeof(double) * p->size);
	  trig0 = (double *) bench_malloc(sizeof(double) * 2 * p->n[0]);
	  if (p->sign < 0)
	       initrc(p->n[0], (double*)p->in, trig0, work);
	  else
	       initcr(p->n[0], (double*)p->in, trig0, work);
     }
     else {
	  double *x = (double*)p->in, *y = x + p->size;
	  work = (double*) bench_malloc(sizeof(double) * p->size * 2);

	  switch (p->rank) {
	      case 1:
		   trig0 =
			(double *) bench_malloc(sizeof(double) * 2 * p->n[0]);
		   init(p->n[0], x, y, trig0, work);
		   break;
	      case 2:
		   trig0 = 
			(double *) bench_malloc(sizeof(double) * 2 * p->n[0]);
		   trig1 =
			(double *) bench_malloc(sizeof(double) * 2 * p->n[1]);
		   init2(p->n[0], p->n[1], x, y, trig0, trig1, work);
		   break;
	      case 3:
		   trig0 =
			(double *) bench_malloc(sizeof(double) * 2 * p->n[0]);
		   trig1 = 
			(double *) bench_malloc(sizeof(double) * 2 * p->n[1]);
		   trig2 = 
			(double *) bench_malloc(sizeof(double) * 2 * p->n[2]);
		   init3(p->n[0], p->n[1], p->n[2], x, y,
			 trig0, trig1, trig2, work);
		   break;
	  }
     }
}

void doit(int iter, struct problem *p)
{
     int i;

     if (p->kind == PROBLEM_COMPLEX) {
	  double *x, *y;
	  
	  if (p->sign < 0) {
	       x = (double *) p->in;
	       y = x + p->size;
	  }
	  else {  /* swap real/imag for inverse transform */
	       y = (double *) p->in;
	       x = y + p->size;
	  }
	  
	  switch (p->rank) {
	      case 1:
	      {
		   unsigned int n0 = p->n[0];
		   for (i = 0; i < iter; ++i) {
			fft(n0, x, y, trig0, work);
		   }
		   break;
	      }
	      case 2:
	      {
		   unsigned int n0 = p->n[0], n1 = p->n[1];
		   for (i = 0; i < iter; ++i) {
			fft2(n0, n1, x, y, trig0, trig1, work);
		   }
		   break;
	      }
	      case 3:
	      {
		   unsigned int n0 = p->n[0], n1 = p->n[1], n2 = p->n[2];
		   for (i = 0; i < iter; ++i) {
			fft3(n0, n1, n2, x, y, trig0, trig1, trig2, work);
		   }
		   break;
	      }
	  }
     }
     else /* PROBLEM_REAL */ {
	  double *x = (double *) p->in;
	  unsigned int n0 = p->n[0];

	  if (p->sign < 0)
	       for (i = 0; i < iter; ++i) {
		    fftrc(n0, x, trig0, work);
	       }
	  else
	       for (i = 0; i < iter; ++i) {
		    fftcr(n0, x, trig0, work);
	       }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(trig1);
     bench_free(trig0);
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

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_1d_halfcomplex(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_halfcomplex(p, in, -1.0);

     /* Need to conjugate input for backwards (halfcomplex->real) transform,
	since both forwards and backwards NAG transforms have the same sign. */
     if (p->sign > 0)
	  conj_hc(p->n[0], (double *) p->in);
}
