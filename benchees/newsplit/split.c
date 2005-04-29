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

#include "split.h"

/* keep a count of the non-trivial real additions and multiplications */
int adds = 0;
int muls = 0;

int twids = 0; /* count of non-trivial (real) twiddle loads */

Ct omega(int n, int m)
{
     Rt theta = (-TWOPI / n) * m;
     return (cos(theta) + I * sin(theta));
}

/* Really, this should be inlined so that when it is multiplied by
   something the * 1.0 multiplication can be eliminated. */
Ct omega_tan(int n, int m)
{
     Rt theta = (-TWOPI / n) * m;
     return (1.0 + I * tan(theta));
}

Rt wcos(int n, int m)
{
     return cos((TWOPI / n) * m);
}
