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

/* $Id: verify.c,v 1.13 2001-07-25 16:31:50 athena Exp $ */

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

static double cerror(bench_complex *A, bench_complex *B, unsigned int n)
{
     /* compute the relative error */
     double error = 0.0;
     double tol = tolerance();
     unsigned int i;

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
static void arand(bench_complex *A, unsigned int n)
{
     unsigned int i;

     for (i = 0; i < n; ++i) {
	  c_re(A[i]) = bench_drand();
	  c_im(A[i]) = bench_drand();
     }
}

/* make array real */
static void mkreal(bench_complex *A, unsigned int n)
{
     unsigned int i;

     for (i = 0; i < n; ++i) {
	  c_im(A[i]) = 0.0;
     }
}

/* make array hermitian */
static void mkhermitian(bench_complex *A, unsigned int n)
{
     unsigned int i;

     c_im(A[0]) = 0.0;

     for (i = 1; 2 * i < n; ++i) {
	  c_re(A[n - i]) = c_re(A[i]);
	  c_im(A[n - i]) = -c_im(A[i]);
     }
     
     if (2 * i == n) {
	  c_im(A[i]) = 0.0;
     }
}


/* C = A - B */
static void asub(bench_complex *C, bench_complex *A, bench_complex *B, unsigned int n)
{
     unsigned int i;

     for (i = 0; i < n; ++i) {
	  c_re(C[i]) = c_re(A[i]) - c_re(B[i]);
	  c_im(C[i]) = c_im(A[i]) - c_im(B[i]);
     }
}

/* B = rotate left A */
static void arol(bench_complex *B, bench_complex *A,
		 unsigned int n, unsigned int n_before, unsigned int n_after)
{
     unsigned int i, ib, ia;

     for (ib = 0; ib < n_before; ++ib) {
	  for (i = 0; i < n - 1; ++i)
	       for (ia = 0; ia < n_after; ++ia)
		    B[(ib * n + i) * n_after + ia] =
			 A[(ib * n + i + 1) * n_after + ia];

	  for (ia = 0; ia < n_after; ++ia)
	       B[(ib * n + n - 1) * n_after + ia] = A[ib * n * n_after + ia];
     }
}

static void aphase_shift(bench_complex *B, bench_complex *A,
			 unsigned int n, 
			 unsigned int n_before, unsigned int n_after,
			 bench_real sign)
{
     unsigned int j, jb, ja;
     const double k2pi = 6.2831853071795864769252867665590057683943388;
     double twopin;
     twopin = k2pi / n;

     for (jb = 0; jb < n_before; ++jb)
	  for (j = 0; j < n; ++j) {
	       double s = sign * sin(j * twopin);
	       double c = cos(j * twopin);

	       for (ja = 0; ja < n_after; ++ja) {
		    unsigned int index = (jb * n + j) * n_after + ja;
		    c_re(B[index]) = c_re(A[index]) * c - c_im(A[index]) * s;
		    c_im(B[index]) = c_re(A[index]) * s + c_im(A[index]) * c;
	       }
	  }
}

static void acmp(bench_complex *A, bench_complex *B, unsigned int n, 
		 const char *test)
{
     double d = cerror(A, B, n);
     if (d > tolerance()) {
	  printf("Found relative error %e (%s)\n", d, test);
	  exit(EXIT_FAILURE);
     }
}

static void do_fft(struct problem *p, bench_complex *in, bench_complex *out)
{
     problem_ccopy_from(p, in);
     doit(1, p);
     problem_ccopy_to(p, out);
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
		   unsigned int rounds)
{
     unsigned int N = p->size;
     unsigned int i;

     /* test 1: check linearity */
     for (i = 0; i < rounds; ++i) {
	  bench_complex alpha, beta;
	  c_re(alpha) = bench_drand();
	  c_re(beta) = bench_drand();

	  if (p->kind == PROBLEM_COMPLEX) {
	       c_im(alpha) = bench_drand();
	       c_im(beta) = bench_drand();
	  } else {
	       c_im(alpha) = c_im(beta) = 0;
	  }
	  arand(inA, N);
	  arand(inB, N);
	  do_fft(p, inA, outA);
	  do_fft(p, inB, outB);

	  cascale(outA, N, alpha);
	  cascale(outB, N, beta);
	  caadd(tmp, outA, outB, N);
	  cascale(inA, N, alpha);
	  cascale(inB, N, beta);
	  caadd(inC, inA, inB, N);
	  do_fft(p, inC, outC);

	  acmp(outC, tmp, N, "linearity");
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
		    unsigned int rounds)
{
     unsigned int n = p->size;
     const bench_complex one = {1.0, 0.0};
     const bench_complex zero = {0.0, 0.0};
     unsigned int i;

     /* test 2: check that the unit impulse is transformed properly */
     caset(inA, n, zero);
     inA[0] = one;
     caset(outA, n, one);

     /* a simple test first, to help with debugging: */
     do_fft(p, inA, outB);
     acmp(outB, outA, n, "impulse response");

     for (i = 0; i < rounds; ++i) {
	  arand(inB, n);
	  asub(inC, inA, inB, n);
	  do_fft(p, inB, outB);
	  do_fft(p, inC, outC);
	  caadd(tmp, outB, outC, n);
	  acmp(tmp, outA, n, "impulse response");
     }
}

enum { TIME_SHIFT, FREQ_SHIFT };

static void tf_shift(struct problem *p,
		     bench_complex *inA,
		     bench_complex *inB,
		     bench_complex *outA,
		     bench_complex *outB,
		     bench_complex *tmp,
		     unsigned int rounds,
		     int which_shift)
{
     double sign;
     unsigned int n, n_before, n_after, dim;
     unsigned int i;

     n = p->size;
     sign = p->sign;

     /* test 3: check the time-shift property */
     /* the paper performs more tests, but this code should be fine too */

     n_before = 1;
     n_after = n;
     for (dim = 0; dim < p->rank; ++dim) {
	  unsigned int n_cur = p->n[dim];

	  n_after /= n_cur;

	  for (i = 0; i < rounds; ++i) {
	       arand(inA, n);

	       if (which_shift == TIME_SHIFT) {
		    if (p->kind == PROBLEM_REAL) 
			 mkreal(inA, n);

		    arol(inB, inA, n_cur, n_before, n_after);
		    do_fft(p, inA, outA);
		    do_fft(p, inB, outB);
		    aphase_shift(tmp, outB, n_cur, n_before, n_after, sign);
		    acmp(tmp, outA, n, "time shift");
	       } else {
		    if (p->kind == PROBLEM_REAL) 
			 mkhermitian(inA, n);

		    aphase_shift(inB, inA, n_cur, n_before, n_after, -sign);
		    do_fft(p, inA, outA);
		    do_fft(p, inB, outB);
		    arol(tmp, outB, n_cur, n_before, n_after);
		    acmp(tmp, outA, n, "freq shift");
	       }
	  }

	  n_before *= n_cur;
     }
}

static void do_verify(struct problem *p, unsigned int rounds)
{
     bench_complex *inA, *inB, *inC, *outA, *outB, *outC, *tmp;
     unsigned int n = p->size;

     if (rounds == 0)
	  rounds = 20;  /* default value */

     inA = (bench_complex *) bench_malloc(n * sizeof(bench_complex));
     inB = (bench_complex *) bench_malloc(n * sizeof(bench_complex));
     inC = (bench_complex *) bench_malloc(n * sizeof(bench_complex));
     outA = (bench_complex *) bench_malloc(n * sizeof(bench_complex));
     outB = (bench_complex *) bench_malloc(n * sizeof(bench_complex));
     outC = (bench_complex *) bench_malloc(n * sizeof(bench_complex));
     tmp = (bench_complex *) bench_malloc(n * sizeof(bench_complex));

     linear(p, inA, inB, inC, outA, outB, outC, tmp, rounds);
     impulse(p, inA, inB, inC, outA, outB, outC, tmp, rounds);

     if (p->kind == PROBLEM_COMPLEX || p->sign == -1)
	  tf_shift(p, inA, inB, outA, outB, tmp, rounds, TIME_SHIFT);
     if (p->kind == PROBLEM_COMPLEX || p->sign == 1)
	  tf_shift(p, inA, inB, outA, outB, tmp, rounds, FREQ_SHIFT);

     bench_free(tmp);
     bench_free(outC);
     bench_free(outB);
     bench_free(outA);
     bench_free(inC);
     bench_free(inB);
     bench_free(inA);
}

void verify(const char *param, int rounds)
{
     struct problem *p;
     p = problem_parse(param);
     BENCH_ASSERT(can_do(p));
     problem_zero(p);
     setup(p);
     do_verify(p, rounds);
     problem_destroy(p);
     if (verbose)
	  ovtpvt("#t\n");
}
