/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

#include "doit-doc.h"

static double *t, *w, ***a;
static int *ip;

#ifdef FORTRAN
   extern void F77_FUNC(cdft3d,CDFT3D)(unsigned int *n1max,
				       unsigned int *n2max,
				       unsigned int *n1, unsigned int *n2,
				       unsigned int *n3,
				       int *isgn, double *a, double *t, 
				       int *ip, double *w);
   extern void F77_FUNC(rdft3d,RDFT3D)(unsigned int *n1max,
				       unsigned int *n2max,
				       unsigned int *n1, unsigned int *n2,
				       unsigned int *n3,
				       int *isgn, double *a, double *t, 
				       int *ip, double *w);

#  define CDFT3D(n1,n2,n3, isgn, a, t,ip,w) F77_FUNC(cdft3d,CDFT3D) \
      (&n3, &n2, &n3, &n2, &n1, &isgn, a[0][0], t,ip,w)
#  define RDFT3D(n1,n2,n3, isgn, a, t,ip,w) F77_FUNC(rdft3d,RDFT3D) \
      (&n3, &n2, &n3, &n2, &n1, &isgn, a[0][0], t,ip,w)

#else
   extern void cdft3d(int n1, int n2, int n3, int isgn, double ***a,
		      double *t, int *ip, double *w);
   extern void rdft3d(int n1, int n2, int n3, int isgn, double ***a,
		      double *t, int *ip, double *w);
#  define CDFT3D(n1,n2,n3, isgn, a, t,ip,w) cdft3d(n1,n2,n3, isgn, a, t,ip,w)
#  define RDFT3D(n1,n2,n3, isgn, a, t,ip,w) rdft3d(n1,n2,n3, isgn, a, t,ip,w)
#endif

int can_do(struct problem *p)
{
     return (DOUBLE_PRECISION &&
	     p->rank == 3 &&
	     (p->n[2] >= 2 || (p->kind == PROBLEM_COMPLEX && p->n[2] >= 1)) &&
	     problem_power_of_two(p, 1));
}

#define IDX(i,j,k) (((i) * n2 + (j)) * n3 + (k))
#define C(i,n) (((n)-(i))%(n))

void copy_h2c(struct problem *p, bench_complex *out)
{
     unsigned int k1, k2, k3, n1, n2, n3;

     n1 = p->n[0];
     n2 = p->n[1];
     n3 = p->n[2];

#define R(i,j,k) c_re(out[IDX(i,j,k)])
#define I(i,j,k) c_im(out[IDX(i,j,k)])

     /* This is ridiculous. */

     for (k1 = 0; k1 < n1; ++k1)
	  for (k2 = 0; k2 < n2; ++k2)
	       for (k3 = 1; k3 < n3/2; ++k3) {
		    R(k1,k2,k3) = a[k1][k2][2*k3];
		    R(C(k1,n1),C(k2,n2),n3-k3) = a[k1][k2][2*k3];
		    I(k1,k2,k3) = -a[k1][k2][2*k3+1];
		    I(C(k1,n1),C(k2,n2),n3-k3) = a[k1][k2][2*k3+1];
	       }

     for (k1 = 0; k1 < n1; ++k1)
	  for (k2 = 1; k2 < n2/2; ++k2) {
	       R(k1,k2,0) = a[k1][k2][0];
	       R(C(k1,n1),n2-k2,0) = a[k1][k2][0];
	       I(k1,k2,0) = -a[k1][k2][1];
	       I(C(k1,n1),n2-k2,0) = a[k1][k2][1];

	       /* Note that Ooura's 2001/11/22 version incorrectly
		  documents the following elements (the first index is
		  swapped). */
	       R(C(k1,n1),k2,n3/2) = a[k1][n2-k2][1];
	       R(k1,n2-k2,n3/2) = a[k1][n2-k2][1];
	       I(C(k1,n1),k2,n3/2) = a[k1][n2-k2][0];
	       I(k1,n2-k2,n3/2) = -a[k1][n2-k2][0];
	  }

     for (k1 = 1; k1 < n1/2; ++k1) {
	  R(k1,0,0) = a[k1][0][0];
	  R(n1-k1,0,0) = a[k1][0][0];
	  I(k1,0,0) = -a[k1][0][1];
	  I(n1-k1,0,0) = a[k1][0][1];
	  R(k1,n2/2,0) = a[k1][n2/2][0];
	  R(n1-k1,n2/2,0) = a[k1][n2/2][0];
	  I(k1,n2/2,0) = -a[k1][n2/2][1];
	  I(n1-k1,n2/2,0) = a[k1][n2/2][1];
	  R(k1,0,n3/2) = a[n1-k1][0][1];
	  R(n1-k1,0,n3/2) = a[n1-k1][0][1];
	  I(k1,0,n3/2) = a[n1-k1][0][0];
	  I(n1-k1,0,n3/2) = -a[n1-k1][0][0];
	  R(k1,n2/2,n3/2) = a[n1-k1][n2/2][1];
	  R(n1-k1,n2/2,n3/2) = a[n1-k1][n2/2][1];
	  I(k1,n2/2,n3/2) = a[n1-k1][n2/2][0];
	  I(n1-k1,n2/2,n3/2) = -a[n1-k1][n2/2][0];
     }

     R(0,0,0) = a[0][0][0];
     R(0,0,n3/2) = a[0][0][1];
     R(0,n2/2,0) = a[0][n2/2][0];
     R(0,n2/2,n3/2) = a[0][n2/2][1];
     R(n1/2,0,0) = a[n1/2][0][0];
     R(n1/2,0,n3/2) = a[n1/2][0][1];
     R(n1/2,n2/2,0) = a[n1/2][n2/2][0];
     R(n1/2,n2/2,n3/2) = a[n1/2][n2/2][1];

     I(0,0,0) = 0;
     I(0,0,n3/2) = 0;
     I(0,n2/2,0) = 0;
     I(0,n2/2,n3/2) = 0;
     I(n1/2,0,0) = 0;
     I(n1/2,0,n3/2) = 0;
     I(n1/2,n2/2,0) = 0;
     I(n1/2,n2/2,n3/2) = 0;
     
#undef R
#undef I
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     unsigned int k1, k2, k3, n1, n2, n3;

     n1 = p->n[0];
     n2 = p->n[1];
     n3 = p->n[2];

#define R(i,j,k) c_re(in[IDX(i,j,k)])
#define I(i,j,k) c_im(in[IDX(i,j,k)])

     for (k1 = 0; k1 < n1; ++k1)
	  for (k2 = 0; k2 < n2; ++k2)
	       for (k3 = 1; k3 < n3/2; ++k3) {
		    a[k1][k2][2*k3] = R(k1,k2,k3);
		    a[k1][k2][2*k3+1] = I(C(k1,n1),C(k2,n2),n3-k3);
	       }

     for (k1 = 0; k1 < n1; ++k1)
	  for (k2 = 1; k2 < n2/2; ++k2) {
	       a[k1][k2][0] = R(k1,k2,0);
	       a[k1][k2][1] = I(C(k1,n1),n2-k2,0); 
	       a[k1][n2-k2][1] = R(C(k1,n1),k2,n3/2);
	       a[k1][n2-k2][0] = I(C(k1,n1),k2,n3/2);
	  }

     for (k1 = 1; k1 < n1/2; ++k1) {
	  a[k1][0][0] = R(k1,0,0);
	  a[k1][0][1] = I(n1-k1,0,0);
	  a[k1][n2/2][0] = R(k1,n2/2,0);
	  a[k1][n2/2][1] = I(n1-k1,n2/2,0);
	  a[n1-k1][0][1] = R(k1,0,n3/2);
	  a[n1-k1][0][0] = I(k1,0,n3/2);
	  a[n1-k1][n2/2][1] = R(k1,n2/2,n3/2);
	  a[n1-k1][n2/2][0] = I(k1,n2/2,n3/2);
     }

     a[0][0][0] = R(0,0,0);
     a[0][0][1] = R(0,0,n3/2);
     a[0][n2/2][0] = R(0,n2/2,0);
     a[0][n2/2][1] = R(0,n2/2,n3/2);
     a[n1/2][0][0] = R(n1/2,0,0);
     a[n1/2][0][1] = R(n1/2,0,n3/2);
     a[n1/2][n2/2][0] = R(n1/2,n2/2,0);
     a[n1/2][n2/2][1] = R(n1/2,n2/2,n3/2);
     
#undef R
#undef I

     ascale(p->in, p->size, 2.0);
}

#undef IDX
#undef C

#define MAX2(a,b) ((a) > (b) ? (a) : (b))
#define MAX3(a,b,c) MAX2(a,MAX2(b,c))

void setup(struct problem *p)
{
     int sip;
     unsigned int i, j, n1, n2, n3, n;
     int sign = 1;
     double *x = (double *) p->in;
 
     BENCH_ASSERT(can_do(p));
     n1 = p->n[0];
     n2 = p->n[1];
     n3 = p->n[2];

     n = p->kind == PROBLEM_COMPLEX ? MAX3(n1,n2,n3) : MAX3(n1,n2,n3/2);

     t = bench_malloc((8 * MAX2(n1,n2)) * sizeof(double));
 
     w = bench_malloc((n/2 + (p->kind == PROBLEM_COMPLEX ? 0 : n3/4))
		       * sizeof(double));

     sip = 2+(1<<(int)(log(n+0.5)/log(2.0))/2);
     ip = bench_malloc(sip * sizeof(int));

     if (p->kind == PROBLEM_COMPLEX)
	  n3 *= 2;
     a = (double ***) bench_malloc(n1 * sizeof(double **));
     for (i = 0; i < n1; ++i) {
	  a[i] = (double **) bench_malloc(n2 * sizeof(double *));
	  for (j = 0; j < n2; ++j)
	       a[i][j] = x + (i * n2 + j) * n3;
     }

     ip[0] = 0;
     if (p->kind == PROBLEM_COMPLEX) {
	  CDFT3D(n1, n2, n3, sign, a, t,ip,w);
     } else {
	  sign = -sign;  /* +1 ==> real->complex */
	  RDFT3D(n1, n2, n3, sign, a, t,ip,w);
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     unsigned int n1 = p->n[0];
     unsigned int n2 = p->n[1];
     unsigned int n3 = p->n[2];
     int sign = p->sign;
     double ***A = a;
     double *T = t;
     double *W = w;
     int *IP = ip;

     if (p->kind == PROBLEM_COMPLEX) {
	  n3 *= 2;  /* grrr */
	  for (i = 0; i < iter; ++i) {
	       CDFT3D(n1, n2, n3, sign, A, T,IP,W);
	  }
     } else {
	  sign = -sign; /* +1 ==> real->complex */
	  for (i = 0; i < iter; ++i) {
	       RDFT3D(n1, n2, n3, sign, A, T,IP,W);
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(a);
     bench_free(t);
     bench_free(w);
     bench_free(ip);
}
