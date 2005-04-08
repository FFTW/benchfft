/* Copyright (c) 2005 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef SPLIT_H
#define SPLIT_H

#include "bench-user.h"

#include <math.h>
#include <stdlib.h>

/* require C99 complex-number support (e.g. use GNU gcc and glibc) */
#include <complex.h>

#if defined(BENCHFFT_SINGLE)
typedef complex float C;
#elif defined(BENCHFFT_LDOUBLE)
typedef complex long double C;
#else
typedef complex double C;
#endif
typedef bench_real R;

#ifdef BENCHFFT_LDOUBLE
typedef long double trigreal;
#  define cos cosl
#  define sin sinl
#  define tan tanl
#  define KTRIG(x) (x##L)
#else
typedef double trigreal;
#  define KTRIG(x) (x)
#endif

#define TWOPI KTRIG(6.2831853071795864769252867665590057683943388)

#define COMPLEX_MUL_3x3 0 /* 1 to count ops assuming 3+3 complex multiply */

/* number of additions and multiplications for a generic complex multiply: */
#if COMPLEX_MUL_3x3
#  define COMPLEX_MUL_ADDS 3
#  define COMPLEX_MUL_MULS 3
#else
#  define COMPLEX_MUL_ADDS 2
#  define COMPLEX_MUL_MULS 4
#endif

/* keep a count of the non-trivial real additions and multiplications */
extern int adds;
extern int muls;

C omega(int n, int m);
C omega_tan(int n, int m);
R wcos(int n, int m);

void oldsplit(int n, C *in0, C *in1, int is, C *out, int os);

R scale(int n, int m);
void nsplit(int n, C *in0, C *in1, int is, C *out, int os);
void nsplitds(int n, C *in0, C *in1, int is, C *out, int os);
void nsplitds2(int n, C *in0, C *in1, int is, C *out, int os);
void nsplitds4(int n, C *in0, C *in1, int is, C *out, int os);

#endif /* SPLIT_H */
