/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

#include "doit-doc.h"

static double *t, *w, **a;
static int *ip;

#ifdef FORTRAN
   extern void F77_FUNC(cdft2d,CDFT2D)(unsigned int *n1max,
				       unsigned int *n1, unsigned int *n2,
				       int *isgn, double *a,
#ifndef NO_T
				       double *t, 
#endif
				       int *ip, double *w);
   extern void F77_FUNC(rdft2d,RDFT2D)(unsigned int *n1max,
				       unsigned int *n1, unsigned int *n2,
				       int *isgn, double *a,
#ifndef NO_T
				       double *t, 
#endif
				       int *ip, double *w);

#  ifndef NO_T
#    define CDFT2D(n1, n2, isgn, a, t, ip, w) F77_FUNC(cdft2d,CDFT2D) \
        (&n2, &n2, &n1, &isgn, a[0], t, ip, w)
#    define RDFT2D(n1, n2, isgn, a, t, ip, w) F77_FUNC(rdft2d,RDFT2D) \
        (&n2, &n2, &n1, &isgn, a[0], t, ip, w)
#  else
#    define CDFT2D(n1, n2, isgn, a, t, ip, w) F77_FUNC(cdft2d,CDFT2D) \
        (&n2, &n2, &n1, &isgn, a[0], ip, w)
#    define RDFT2D(n1, n2, isgn, a, t, ip, w) F77_FUNC(rdft2d,RDFT2D) \
        (&n2, &n2, &n1, &isgn, a[0], ip, w)
#  endif

#else
   extern void cdft2d(int n1, int n2, int isgn, double **a,
#  ifndef NO_T
		      double *t,
#  endif
		      int *ip, double *w);
   extern void rdft2d(int n1, int n2, int isgn, double **a,
#  ifndef NO_T
		      double *t,
#  endif
		      int *ip, double *w);
#  ifndef NO_T
#    define CDFT2D(n1, n2, isgn, a, t, ip, w) cdft2d(n1, n2, isgn, a, t, ip, w)
#    define RDFT2D(n1, n2, isgn, a, t, ip, w) rdft2d(n1, n2, isgn, a, t, ip, w)
#  else
#    define CDFT2D(n1, n2, isgn, a, t, ip, w) cdft2d(n1, n2, isgn, a, ip, w)
#    define RDFT2D(n1, n2, isgn, a, t, ip, w) rdft2d(n1, n2, isgn, a, ip, w)
#  endif
#endif

int can_do(struct problem *p)
{
     return (DOUBLE_PRECISION &&
	     p->rank == 2 &&
	     (p->n[1] >= 2 || (p->kind == PROBLEM_COMPLEX && p->n[1] >= 1)) &&
	     problem_power_of_two(p, 1));
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     unsigned int k1, k2, n1, n2;

     n1 = p->n[0];
     n2 = p->n[1];

#define IDX(i,j) ((i) * n2 + (j))

     for (k1 = 1; k1 < n1; ++k1)
	  for (k2 = 1; k2 < n2/2; ++k2) {
	       c_re(out[IDX(k1,k2)]) = a[k1][2*k2];
	       c_re(out[IDX(n1-k1,n2-k2)]) = a[k1][2*k2];
	       c_im(out[IDX(k1,k2)]) = -a[k1][2*k2 + 1];
	       c_im(out[IDX(n1-k1,n2-k2)]) = a[k1][2*k2 + 1];
	  }

     for (k2 = 1; k2 < n2/2; ++k2) {
	  c_re(out[IDX(0,k2)]) = a[0][2*k2];
	  c_re(out[IDX(0,n2-k2)]) = a[0][2*k2];
	  c_im(out[IDX(0,k2)]) = -a[0][2*k2+1];
	  c_im(out[IDX(0,n2-k2)]) = a[0][2*k2+1];
     }

     for (k1 = 1; k1 < n1/2; ++k1) {
	  c_re(out[IDX(k1,0)]) = a[k1][0];
	  c_re(out[IDX(n1-k1,0)]) = a[k1][0];
	  c_im(out[IDX(k1,0)]) = -a[k1][1];
	  c_im(out[IDX(n1-k1,0)]) = a[k1][1];
	  c_re(out[IDX(k1,n2/2)]) = a[n1-k1][1];
	  c_re(out[IDX(n1-k1,n2/2)]) = a[n1-k1][1];
	  c_im(out[IDX(k1,n2/2)]) = a[n1-k1][0];
	  c_im(out[IDX(n1-k1,n2/2)]) = -a[n1-k1][0];
     }

     c_re(out[IDX(0,0)]) = a[0][0];
     c_re(out[IDX(0,n2/2)]) = a[0][1];
     c_re(out[IDX(n1/2,0)]) = a[n1/2][0];
     c_re(out[IDX(n1/2,n2/2)]) = a[n1/2][1];
     c_im(out[IDX(0,0)]) = 0;
     c_im(out[IDX(0,n2/2)]) = 0;
     c_im(out[IDX(n1/2,0)]) = 0;
     c_im(out[IDX(n1/2,n2/2)]) = 0;

#undef IDX
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     unsigned int k1, k2, n1, n2;

     n1 = p->n[0];
     n2 = p->n[1];

#define IDX(i,j) ((i) * n2 + (j))

     for (k1 = 1; k1 < n1; ++k1)
	  for (k2 = 1; k2 < n2/2; ++k2) {
	       a[k1][2*k2] = c_re(in[IDX(k1,k2)]);
	       a[k1][2*k2 + 1] = c_im(in[IDX(n1-k1,n2-k2)]);
	  }

     for (k2 = 1; k2 < n2/2; ++k2) {
	  a[0][2*k2] = c_re(in[IDX(0,k2)]);
	  a[0][2*k2+1] = c_im(in[IDX(0,n2-k2)]);
     }

     for (k1 = 1; k1 < n1/2; ++k1) {
	  a[k1][0] = c_re(in[IDX(k1,0)]);
	  a[k1][1] = c_im(in[IDX(n1-k1,0)]);
	  a[n1-k1][1] = c_re(in[IDX(k1,n2/2)]);
	  a[n1-k1][0] = c_im(in[IDX(k1,n2/2)]);
     }

     a[0][0] = c_re(in[IDX(0,0)]);
     a[0][1] = c_re(in[IDX(0,n2/2)]);
     a[n1/2][0] = c_re(in[IDX(n1/2,0)]);
     a[n1/2][1] = c_re(in[IDX(n1/2,n2/2)]);

#undef IDX

     ascale(p->in, p->size, 2.0);
}

void setup(struct problem *p)
{
     int sip;
     unsigned int i, n1, n2, n = 0;
     int sign = 1;
     double *x = (double *) p->in;
 
     BENCH_ASSERT(can_do(p));
     n1 = p->n[0];
     n2 = p->n[1];
     if (p->kind == PROBLEM_COMPLEX)
	  n = (n1 > n2) ? n1 : n2;
     else
	  n = (n1 > n2/2) ? n1 : n2/2;

#ifndef NO_T
     t = bench_malloc((8 * n1) * sizeof(double));
#endif
 
     w = bench_malloc((n/2 + (p->kind == PROBLEM_COMPLEX ? 0 : n2/4))
		       * sizeof(double));

     sip = 2+(1<<(int)(log(n+0.5)/log(2.0))/2);
     ip = bench_malloc(sip * sizeof(int));

     if (p->kind == PROBLEM_COMPLEX)
	  n2 *= 2;
     a = (double **) bench_malloc(n1 * sizeof(double *));
     for (i = 0; i < n1; ++i)
	  a[i] = x + i * n2;

     ip[0] = 0;
     if (p->kind == PROBLEM_COMPLEX) {
	  CDFT2D(n1, n2, sign, a, t, ip, w);
     } else {
	  sign = -sign;  /* +1 ==> real->complex */
	  RDFT2D(n1, n2, sign, a, t, ip, w);
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     unsigned int n1 = p->n[0];
     unsigned int n2 = p->n[1];
     int sign = p->sign;
     double **A = a;
     double *T = t;
     double *W = w;
     int *IP = ip;

     if (p->kind == PROBLEM_COMPLEX) {
	  n2 *= 2;  /* grrr */
	  for (i = 0; i < iter; ++i) {
	       CDFT2D(n1, n2, sign, A, T, IP, W);
	  }
     } else {
	  sign = -sign; /* +1 ==> real->complex */
	  for (i = 0; i < iter; ++i) {
	       RDFT2D(n1, n2, sign, A, T, IP, W);
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(a);
#ifndef NO_T
     bench_free(t);
#endif
     bench_free(w);
     bench_free(ip);
}
