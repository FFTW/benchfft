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
BENCH_DOC("name", NAME)
BENCH_DOC("package", "GNU Scientific Library (GSL)")
BENCH_DOC("author", "Brian Gough")
BENCH_DOC("year", "2002")
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


#define CAT0(a, b) a ## b
#define CAT(a, b) CAT0(a, b)

#ifdef BENCHFFT_SINGLE
#include <gsl/gsl_fft_complex_float.h>
#include <gsl/gsl_fft_real_float.h>
#include <gsl/gsl_fft_halfcomplex_float.h>
#define C(x) CAT(gsl_fft_complex_float_, P(x))
#define R(x) CAT(gsl_fft_real_float_radix2_, x)
#define H(x) CAT(gsl_fft_halfcomplex_float_radix2_, x)
#else
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#define C(x) CAT(gsl_fft_complex_, P(x))
#define R(x) CAT(gsl_fft_real_radix2_, x)
#define H(x) CAT(gsl_fft_halfcomplex_radix2_, x)
#endif


int can_do(struct problem *p)
{
     return (p->rank == 1 &&
#ifdef DO_COMPLEX
	     p->kind == PROBLEM_COMPLEX
#else
	     p->kind == PROBLEM_REAL
#endif
	     && problem_power_of_two(p, 1));
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_1d_halfcomplex(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_halfcomplex(p, in, -1.0);
}

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
     /* nothing to do */
}

void doit(int iter, struct problem *p)
{
     int i;
     int n = p->n[0];
     void *in = p->in;

     if (p->kind == PROBLEM_COMPLEX) {
#ifdef DO_COMPLEX
	  if (p->sign == -1) {
	    for (i = 0; i < iter; ++i) 
	      C(forward)(in, 1, n);
	  } else {
	    for (i = 0; i < iter; ++i) 
	      C(backward)(in, 1, n);
	  }
#endif
     } else {
#ifdef DO_REAL
	  if (p->sign == -1) {
	    for (i = 0; i < iter; ++i) 
	      R(transform)(in, 1, n);
	  } else {
	    for (i = 0; i < iter; ++i) 
	      H(transform)(in, 1, n);
	  }
#endif
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}

