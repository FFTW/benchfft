/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>
#include <stdio.h>

#ifdef HAVE_SUNPERF_H
#include <sunperf.h>
#endif

#include "common.h"

#ifdef BENCHFFT_SINGLE
#define CFFTC2 cfftc2_
#define SFFTC2 sfftc2_
#define CFFTS2 cffts2_
#else
#define CFFTC2 zfftz2_
#define SFFTC2 dfftz2_
#define CFFTS2 zfftd2_
#endif

int can_do(struct problem *p)
{
     return (p->rank == 2);
}

static void *WORK;
static int lwork;
static void *TRIGS;
static int IFAC[2*128];
static int ldx;
static int ldy;

void setup(struct problem *p)
{
     int iopt, n1, n2, ierr;
     bench_real scale = 1.0;

     BENCH_ASSERT(can_do(p));
     n1 = p->n[1];
     n2 = p->n[0];

     TRIGS = bench_malloc((2 * (n1 + n2)) * sizeof(bench_real));
     iopt = 0;
 
     if (p->kind == PROBLEM_COMPLEX) {
	  ldx = ldy = n1;
	  lwork = 2 * MAX2(n1, 2*n2);
	  WORK = bench_malloc(lwork * sizeof(bench_real));
	  CFFTC2(&iopt, &n1, &n2, &scale, 0, &ldx, 0, &ldy,
		 TRIGS, IFAC, WORK, &lwork, &ierr);
     } else {
	  lwork = MAX2(n1, 2*n2);
	  WORK = bench_malloc(lwork * sizeof(bench_real));
	  if (p->sign == -1) {
	       ldy = n1 / 2 + 1;
	       ldx = 2 * ldy;
	       SFFTC2(&iopt, &n1, &n2, &scale, 0, &ldx, 0, &ldy,
		      TRIGS, IFAC, WORK, &lwork, &ierr);
	  } else {
	       ldx = n1 / 2 + 1;
	       ldy = 2 * ldx;
	       CFFTS2(&iopt, &n1, &n2, &scale, 0, &ldx, 0, &ldy,
		      TRIGS, IFAC, WORK, &lwork, &ierr);
	  }
     }
     BENCH_ASSERT(ierr == 0);
}

void doit(int iter, struct problem *p)
{
     int i;
     int n1 = p->n[1];
     int n2 = p->n[0];
     int iopt = p->sign;
     void *in = p->in;
     void *out = p->out;
     bench_real scale = 1.0;
     int ierr;

     if (p->kind == PROBLEM_COMPLEX) {
	  for (i = 0; i < iter; ++i) 
	       CFFTC2(&iopt, &n1, &n2, &scale, in, &ldx, out, &ldy,
		      TRIGS, IFAC, WORK, &lwork, &ierr);
     } else {
	  if (iopt == -1) {
	       for (i = 0; i < iter; ++i) {
		    SFFTC2(&iopt, &n1, &n2, &scale, in, &ldx, out, &ldy,
			   TRIGS, IFAC, WORK, &lwork, &ierr);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) 
		    CFFTS2(&iopt, &n1, &n2, &scale, in, &ldx, out, &ldy,
			   TRIGS, IFAC, WORK, &lwork, &ierr);
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(WORK);
     bench_free(TRIGS);
}
