/*
 * Copyright (c) 2000 Matteo Frigo
 * Copyright (c) 2000 Steven G. Johnson
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* $Id: verify.c,v 1.2 2001-07-05 16:49:43 athena Exp $ */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "config.h"
#include "bench.h"

/*************************************************
 * complex correctness test
 *************************************************/
static double tolerance(void)
{
     if (sizeof(bench_real) == sizeof(float))
	  return 1.0e-2;
     else
	  return 1.0e-6;
}

#ifndef HAVE_HYPOT
static double hypot(double a, double b)
{
     return sqrt(a * a + b * b);
}
#endif

static double cerror(bench_complex *A, bench_complex *B, int n)
{
     /* compute the relative error */
     double error = 0.0;
     double tol = tolerance();
     int i;

     for (i = 0; i < n; ++i) {
	  double a;
	  double mag;
	  a = hypot(c_re(A[i]) - c_re(B[i]), c_im(A[i]) - c_im(B[i]));
	  mag = 0.5 * (hypot(c_re(A[i]), c_im(A[i])) + 
		       hypot(c_re(B[i]), c_im(B[i])))
	       + tol;

	  a /= mag;
	  if (a > error)
	       error = a;

#ifdef HAVE_ISNAN
	  BENCH_ASSERT(!isnan(a));
#endif
     }
     return error;
}

/* generate random inputs */
static void arand(bench_complex *A, int n)
{
     int i;

     for (i = 0; i < n; ++i) {
	  c_re(A[i]) = bench_drand();
	  c_im(A[i]) = bench_drand();
     }
}

/* C = A + B */
static void aadd(bench_complex *C, bench_complex *A, bench_complex *B, int n)
{
     int i;

     for (i = 0; i < n; ++i) {
	  c_re(C[i]) = c_re(A[i]) + c_re(B[i]);
	  c_im(C[i]) = c_im(A[i]) + c_im(B[i]);
     }
}


/* C = A - B */
static void asub(bench_complex *C, bench_complex *A, bench_complex *B, int n)
{
     int i;

     for (i = 0; i < n; ++i) {
	  c_re(C[i]) = c_re(A[i]) - c_re(B[i]);
	  c_im(C[i]) = c_im(A[i]) - c_im(B[i]);
     }
}

/* B = rotate left A */
static void arol(bench_complex *B, bench_complex *A,
		 int n, int n_before, int n_after)
{
     int i, ib, ia;

     for (ib = 0; ib < n_before; ++ib) {
	  for (i = 0; i < n - 1; ++i)
	       for (ia = 0; ia < n_after; ++ia)
		    B[(ib * n + i) * n_after + ia] =
			 A[(ib * n + i + 1) * n_after + ia];

	  for (ia = 0; ia < n_after; ++ia)
	       B[(ib * n + n - 1) * n_after + ia] = A[ib * n * n_after + ia];
     }
}

/* A = alpha * A  (in place) */
static void ascale(bench_complex *A, bench_complex alpha, int n)
{
     int i;

     for (i = 0; i < n; ++i) {
	  bench_complex a = A[i];
	  c_re(A[i]) = c_re(a) * c_re(alpha) - c_im(a) * c_im(alpha);
	  c_im(A[i]) = c_re(a) * c_im(alpha) + c_im(a) * c_re(alpha);
     }
}

static void acmp(bench_complex *A, bench_complex *B, int n)
{
     double d = cerror(A, B, n);
     if (d > tolerance()) {
	  fflush(stdout);
	  fprintf(stderr, "Found relative error %e\n", d);
	  BENCH_ASSERT(((void)"failure in Ergun's verification procedure\n",0));
     }
}

static void do_fft(struct problem *p, bench_complex *in, bench_complex *out)
{
     cacopy(in, p->p.complex.in, p->size);
     doit(1, p);
     cacopy(p->p.complex.out, out, p->size);
}

/*
 * Implementation of the FFT tester described in
 *
 * Funda Ergün. Testing multivariate linear functions: Overcoming the
 * generator bottleneck. In Proceedings of the Twenty-Seventh Annual
 * ACM Symposium on the Theory of Computing, pages 407-416, Las Vegas,
 * Nevada, 29 May--1 June 1995.
 */

static void linear(struct problem *p,
		   bench_complex *inA,
		   bench_complex *inB,
		   bench_complex *inC,
		   bench_complex *outA,
		   bench_complex *outB,
		   bench_complex *outC,
		   bench_complex *tmp,
		   int rounds)
{
     int N = p->size;
     int i;

     /* test 1: check linearity */
     for (i = 0; i < rounds; ++i) {
	  bench_complex alpha, beta;
	  c_re(alpha) = bench_drand();
	  c_im(alpha) = bench_drand();
	  c_re(beta) = bench_drand();
	  c_im(beta) = bench_drand();
	  arand(inA, N);
	  arand(inB, N);
	  do_fft(p, inA, outA);
	  do_fft(p, inB, outB);

	  ascale(outA, alpha, N);
	  ascale(outB, beta, N);
	  aadd(tmp, outA, outB, N);
	  ascale(inA, alpha, N);
	  ascale(inB, beta, N);
	  aadd(inC, inA, inB, N);
	  do_fft(p, inC, outC);

	  acmp(outC, tmp, N);
     }

}

static void impulse(struct problem *p,
		    bench_complex *inA,
		    bench_complex *inB,
		    bench_complex *inC,
		    bench_complex *outA,
		    bench_complex *outB,
		    bench_complex *outC,
		    bench_complex *tmp,
		    int rounds)
{
     int n = p->size;
     const bench_complex one = {1.0, 0.0};
     const bench_complex zero = {0.0, 0.0};
     int i;

     /* test 2: check that the unit impulse is transformed properly */
     caset(inA, n, zero);
     inA[0] = one;
     caset(outA, n, one);

     /* a simple test first, to help with debugging: */
     do_fft(p, inA, outB);
     acmp(outB, outA, n);

     for (i = 0; i < rounds; ++i) {
	  arand(inB, n);
	  asub(inC, inA, inB, n);
	  do_fft(p, inB, outB);
	  do_fft(p, inC, outC);
	  aadd(tmp, outB, outC, n);
	  acmp(tmp, outA, n);
     }
}


static void do_verify(struct problem *p, int rounds)
{
     bench_complex *inA, *inB, *inC, *outA, *outB, *outC, *tmp;
     int i;
     double sign;
     int n;
     double twopin;
     const double k2pi = 6.2831853071795864769252867665590057683943388;

     /* TODO: real case */
     BENCH_ASSERT(p->kind == PROBLEM_COMPLEX);

     /* TODO: nd */
     BENCH_ASSERT(p->rank == 1);

     if (rounds == 0)
	  rounds = 20;  /* default value */

     n = p->size;
     twopin = k2pi / n;
     sign = p->p.complex.sign;

     inA = (bench_complex *) bench_malloc(n * sizeof(bench_complex));
     inB = (bench_complex *) bench_malloc(n * sizeof(bench_complex));
     inC = (bench_complex *) bench_malloc(n * sizeof(bench_complex));
     outA = (bench_complex *) bench_malloc(n * sizeof(bench_complex));
     outB = (bench_complex *) bench_malloc(n * sizeof(bench_complex));
     outC = (bench_complex *) bench_malloc(n * sizeof(bench_complex));
     tmp = (bench_complex *) bench_malloc(n * sizeof(bench_complex));

     linear(p, inA, inB, inC, outA, outB, outC, tmp, rounds);
     impulse(p, inA, inB, inC, outA, outB, outC, tmp, rounds);

     /* test 3: check the time-shift property */
     /* the paper performs more tests, but this code should be fine too */
     for (i = 0; i < rounds; ++i) {
	  int j;

	  arand(inA, n);
	  arol(inB, inA, n, 1, 1);
	  do_fft(p, inA, outA);
	  do_fft(p, inB, outB);

	  for (j = 0; j < n; ++j) {
	       double s = sign * sin(j * twopin);
	       double c = cos(j * twopin);
	       c_re(tmp[j]) = c_re(outB[j]) * c - c_im(outB[j]) * s;
	       c_im(tmp[j]) = c_re(outB[j]) * s + c_im(outB[j]) * c;
	  }

	  acmp(tmp, outA, n);
     }

     bench_free(tmp);
     bench_free(outC);
     bench_free(outB);
     bench_free(outA);
     bench_free(inC);
     bench_free(inB);
     bench_free(inA);
}

void verify(const char *param)
{
     struct problem *p;
     p = problem_parse(param);
     BENCH_ASSERT(can_do(p));
     setup(p);
     do_verify(p, 10);
     problem_destroy(p);
     printf("ok\n");
}
