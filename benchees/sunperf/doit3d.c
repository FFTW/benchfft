/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>
#include <stdio.h>

#ifdef HAVE_SUNPERF_H
#include <sunperf.h>
#endif

#include "common.h"

#ifdef BENCHFFT_SINGLE
#define CFFTC3 cfftc3_
#define SFFTC3 sfftc3_
#define CFFTS3 cffts3_
#else
#define CFFTC3 zfftz3_
#define SFFTC3 dfftz3_
#define CFFTS3 zfftd3_
#endif

int can_do(struct problem *p)
{
     return (p->rank == 3);
}

static void *WORK;
static int lwork;
static void *TRIGS;
static int IFAC[3*128];
static int ldx1, ldx2;
static int ldy1, ldy2;

void setup(struct problem *p)
{
     int iopt, n1, n2, n3, ierr;
     bench_real scale = 1.0;

     BENCH_ASSERT(can_do(p));
     n1 = p->n[2];
     n2 = p->n[1];
     n3 = p->n[0];

     TRIGS = bench_malloc((2 * (n1 + n2 + n3)) * sizeof(bench_real));
     iopt = 0;
 
     if (p->kind == PROBLEM_COMPLEX) {
	  ldx1 = ldy1 = n1;
	  ldx2 = ldy2 = n2;
	  lwork = 32 * n3 + 2 * MAX3(n1, n2, n3);
	  WORK = bench_malloc(lwork * sizeof(bench_real));
	  CFFTC3(&iopt, &n1, &n2, &n3, &scale, 0, &ldx1, &ldx2,
		 0, &ldy1, &ldy2, TRIGS, IFAC, WORK, &lwork, &ierr);
     } else {
	  lwork = 16 * n3 + 2 * MAX3(n1, 2*n2, 2*n3);
	  WORK = bench_malloc(lwork * sizeof(bench_real));
	  if (p->sign == -1) {
	       ldy1 = n1 / 2 + 1;
	       ldx1 = 2 * ldy1;
	       ldx2 = ldy2 = n2;
	       SFFTC3(&iopt, &n1, &n2, &n3, &scale, 0, &ldx1, &ldx2,
		      0, &ldy1, &ldy2, TRIGS, IFAC, WORK, &lwork, &ierr);
	  } else {
	       ldx1 = n1 / 2 + 1;
	       ldy1 = 2 * ldx1;
	       ldx2 = ldy2 = n2;
	       CFFTS3(&iopt, &n1, &n2, &n3, &scale, 0, &ldx1, &ldx2,
		      0, &ldy1, &ldy2, TRIGS, IFAC, WORK, &lwork, &ierr);
	  }
     }
     BENCH_ASSERT(ierr == 0);
}

void doit(int iter, struct problem *p)
{
     int i;
     int n1 = p->n[2];
     int n2 = p->n[1];
     int n3 = p->n[0];
     int iopt = p->sign;
     void *in = p->in;
     void *out = p->out;
     bench_real scale = 1.0;
     int ierr;

     if (p->kind == PROBLEM_COMPLEX) {
	  for (i = 0; i < iter; ++i) 
	       CFFTC3(&iopt, &n1, &n2, &n3, &scale, in, &ldx1, &ldx2, out,
		      &ldy1, &ldy2, TRIGS, IFAC, WORK, &lwork, &ierr);
     } else {
	  if (iopt == -1) {
	       for (i = 0; i < iter; ++i) {
		    SFFTC3(&iopt, &n1, &n2, &n3, &scale, in, &ldx1, &ldx2, 
			   out, &ldy1, &ldy2, TRIGS, IFAC, WORK, &lwork, 
			   &ierr);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) 
		    CFFTS3(&iopt, &n1, &n2, &n3, &scale, in, &ldx1, &ldx2, 
			   out, &ldy1, &ldy2, TRIGS, IFAC, WORK, &lwork, 
			   &ierr);
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(WORK);
     bench_free(TRIGS);
}
