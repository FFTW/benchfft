/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>
#include <stdio.h>

#ifdef HAVE_SUNPERF_H
#include <sunperf.h>
#endif

#include "common.h"

#ifdef BENCHFFT_SINGLE
#define CFFTC cfftc_
#define SFFTC sfftc_
#define CFFTS cffts_
#else
#define CFFTC zfftz_
#define SFFTC dfftz_
#define CFFTS zfftd_
#endif

int can_do(struct problem *p)
{
     return (p->rank == 1);
}

static void *WORK;
static int lwork;
static void *TRIGS;
static int IFAC[128];

void setup(struct problem *p)
{
     int iopt, n, ierr;
     bench_real scale = 1.0;

     BENCH_ASSERT(can_do(p));
     n = p->n[0];

     lwork = (p->kind == PROBLEM_COMPLEX) ? 2 * n : n;
     WORK = bench_malloc(lwork * sizeof(bench_real));
     TRIGS = bench_malloc((2 * n) * sizeof(bench_real));
     iopt = 0;
 
     if (p->kind == PROBLEM_COMPLEX) {
	  CFFTC(&iopt, &n, &scale, 0, 0, TRIGS, IFAC, WORK, &lwork, &ierr);
     } else {
	  if (p->sign == -1) {
	       SFFTC(&iopt, &n, &scale, 0, 0, TRIGS, IFAC, 
		     WORK, &lwork, &ierr);
	  } else {
	       CFFTS(&iopt, &n, &scale, 0, 0, TRIGS, IFAC, 
		     WORK, &lwork, &ierr);
	  }
     }
     BENCH_ASSERT(ierr == 0);
}

void doit(int iter, struct problem *p)
{
     int i;
     int n = p->n[0];
     int iopt = p->sign;
     void *in = p->in;
     void *out = p->out;
     bench_real scale = 1.0;
     int ierr;

     if (p->kind == PROBLEM_COMPLEX) {
	  for (i = 0; i < iter; ++i) 
	       CFFTC(&iopt, &n, &scale, in, out, TRIGS, IFAC, WORK, 
		     &lwork, &ierr);
     } else {
	  if (iopt == -1) {
	       for (i = 0; i < iter; ++i) {
		    SFFTC(&iopt, &n, &scale, in, out, TRIGS, IFAC, WORK, 
			  &lwork, &ierr);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) 
		    CFFTS(&iopt, &n, &scale, in, out, TRIGS, IFAC, WORK, 
			  &lwork, &ierr);
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(WORK);
     bench_free(TRIGS);
}
