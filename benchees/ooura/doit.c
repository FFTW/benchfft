/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

#include "doit-doc.h"

static double *w;
static int *ip;

#ifdef FORTRAN
   extern void F77_FUNC(cdft,CDFT)(int *n, int *isgn, double *a, 
				   int *ip, double *w);
   extern void F77_FUNC(rdft,RDFT)(int *n, int *isgn, double *a, 
				   int *ip, double *w);

#  define CDFT(n, isgn, a, ip, w) F77_FUNC(cdft,CDFT)(&n, &isgn, a, ip, w)
#  define RDFT(n, isgn, a, ip, w) F77_FUNC(rdft,RDFT)(&n, &isgn, a, ip, w)

#else
   extern void cdft(int n, int isgn, double *a, int *ip, double *w);
   extern void rdft(int n, int isgn, double *a, int *ip, double *w);
#  define CDFT(n, isgn, a, ip, w) cdft(n, isgn, a, ip, w)
#  define RDFT(n, isgn, a, ip, w) rdft(n, isgn, a, ip, w)
#endif

int can_do(struct problem *p)
{
     return (DOUBLE_PRECISION &&
	     p->rank == 1 &&
	     (p->n[0] >= 2 || (p->kind == PROBLEM_COMPLEX && p->n[0] >= 1)) &&
	     problem_power_of_two(p, 1));
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_1d_packed(p, out, 1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_packed(p, in, 1.0);
     ascale(p->in, p->size, 2.0);
}

void setup(struct problem *p)
{
     int sip;
     int n;
     int sign = 1;
 
     BENCH_ASSERT(can_do(p));
     n = p->n[0];
 
     w = bench_malloc((n / 2) * sizeof(double));

     sip = 2+(1<<(int)(log(n+0.5)/log(2.0))/2);
     ip = bench_malloc(sip * sizeof(int));

     ip[0] = 0;
     if (p->kind == PROBLEM_COMPLEX) {
	  n *= 2;
	  CDFT(n, sign, p->in, ip, w);
     } else {
	  sign = -sign;  /* +1 ==> real->complex */
	  RDFT(n, sign, p->in, ip, w);
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     int n = p->n[0];
     int sign = p->sign;
     void *in = p->in;
     double *W = w;
     int *IP = ip;

     if (p->kind == PROBLEM_COMPLEX) {
	  n *= 2;  /* grrr */
	  for (i = 0; i < iter; ++i) {
	       CDFT(n, sign, in, IP, W);
	  }
     } else {
	  sign = -sign; /* +1 ==> real->complex */
	  for (i = 0; i < iter; ++i) {
	       RDFT(n, sign, in, IP, W);
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(w);
     bench_free(ip);
}
