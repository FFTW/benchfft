/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

#ifdef HAVE_GSL_GSL_VERSION_H
#include <gsl/gsl_version.h>
#endif

static const char *mkvers(void)
{
#ifdef HAVE_GSL_GSL_VERSION_H
     return gsl_version;
#else
     return "unknown"
#endif
}

BEGIN_BENCH_DOC
BENCH_DOC("name", "gsl-mixed-radix")
BENCH_DOC("author", "Brian Gough")
BENCH_DOC("year", "2001")
BENCH_DOC("url", "http://sources.redhat.com/gsl")
BENCH_DOC("url-was-valid-on", "Thu Jul 26 07:48:45 EDT 2001")
BENCH_DOC("copyright",
  "Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough\n"
  "\n" 
  "This program is free software; you can redistribute it and/or modify\n"
  "it under the terms of the GNU General Public License as published by\n"
  "the Free Software Foundation; either version 2 of the License, or (at\n"
  "your option) any later version.\n")
BENCH_DOC("language", "C")
BENCH_DOCF("version", mkvers)
END_BENCH_DOC

#ifdef BENCHFFT_SINGLE
#include <gsl/gsl_fft_complex_float.h>
#include <gsl/gsl_fft_real_float.h>
#include <gsl/gsl_fft_halfcomplex_float.h>
#define C(x) gsl_fft_complex_float_##x
#define R(x) gsl_fft_real_float_##x
#define H(x) gsl_fft_halfcomplex_float_##x
#else
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#define C(x) gsl_fft_complex_##x
#define R(x) gsl_fft_real_##x
#define H(x) gsl_fft_halfcomplex_##x
#endif

int can_do(struct problem *p)
{
     return (p->rank == 1 && problem_in_place(p));
}

/* fftpack-style real transforms */
void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_1d_fftpack(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_fftpack(p, in, -1.0);
}

static void *W;

void setup(struct problem *p)
{
     int n;
 
     BENCH_ASSERT(can_do(p));
     n = p->n[0];
 
     if (p->kind == PROBLEM_COMPLEX) {
	  W = C(alloc)(n);
     } else {
	  if (p->sign == -1)
	       W = R(alloc)(n);
	  else
	       W = H(alloc)(n);
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     int n = p->n[0];
     void *in = p->in;
     void *w = W;

     if (p->kind == PROBLEM_COMPLEX) {
	  if (p->sign == -1) {
	    for (i = 0; i < iter; ++i) 
	      C(forward)(in, 1, n, w);
	  } else {
	    for (i = 0; i < iter; ++i) 
	      C(backward)(in, 1, n, w);
	  }
     } else {
	  if (p->sign == -1) {
	    for (i = 0; i < iter; ++i) 
	      R(transform)(in, 1, n, w);
	  } else {
	    for (i = 0; i < iter; ++i) 
	      H(transform)(in, 1, n, w);
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     if (p->kind == PROBLEM_COMPLEX) {
	  C(free)(W);
     } else {
	  if (p->sign == -1)
	       R(free)(W);
	  else
	       H(free)(W);
     }
}

