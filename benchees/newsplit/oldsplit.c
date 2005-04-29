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

/* Ordinary "conjugate-pair" split-radix FFT.
   This code achieves the 1968 operation count of Yavne.

   in0 points to the 0th element of the input, whereas the
   i-th element of in1 (with stride is) is the (i+1)-th element
   of the input.  (Instead of of keeping track of the 0th input
   separately in this way, we could also do the 2nd n/4 sub-transform
   as an inverse (opposite-sign) DFT with a negative stride.) */
void oldsplit(int n, C *in0, C *in1, int is, C *out, int os)
{
     if (n == 1)
	  out[0] = *in0;
     else if (n == 2) {
	  C a, b;
	  a = *in0;
	  b = *in1;
	  out[0] = a + b;
	  out[os] = a - b;

	  adds += 2 * 2;
     }
     else {
	  int i;
	  oldsplit(n/2, in0, in1 + is, is*2, out, os);
	  oldsplit(n/4, in1, in1 + 4*is, is*4, out + os * n/2, os);
	  oldsplit(n/4, in1 + (n-2)*is, in1 + 2*is, is*4, out + os * 3*n/4,os);
	  
	  /* i == 0 case (no multiplies): */
	  {
	       C f0 = out[0*os];
               C f1 = out[(0+n/4)*os];
               C wg = out[(0+n/2)*os];
               C wh = out[(0+3*n/4)*os];
               C wghp = wg + wh;
               C wghm = wg - wh;
               out[0*os] = f0 + wghp;
               out[(0+n/4)*os] = f1 - I * wghm;
               out[(0+n/2)*os] = f0 - wghp;
               out[(0+3*n/4)*os] = f1 + I * wghm;

	       adds += 2 * 6;
	  }

	  /* 1 <= i < n/8 are pairs of related twiddle factors */
	  for (i = 1; i < n/8; ++i) {
	       Ct w = omega(n, i);

	       {
		    C f0 = out[i*os];
		    C f1 = out[(i+n/4)*os];
		    C wg = w * out[(i+n/2)*os];
		    C wh = conj(w) * out[(i+3*n/4)*os];
		    C wghp = wg + wh;
		    C wghm = wg - wh;
		    out[i*os] = f0 + wghp;
		    out[(i+n/4)*os] = f1 - I * wghm;
		    out[(i+n/2)*os] = f0 - wghp;
		    out[(i+3*n/4)*os] = f1 + I * wghm;
	       }

	       {
		    int j = n/4 - i;
		    C f0 = out[j*os];
		    C f1 = out[(j+n/4)*os];
		    C wg = -I * conj(w) * out[(j+n/2)*os];
		    C wh = I * w * out[(j+3*n/4)*os];
		    C wghp = wg + wh;
		    C wghm = wg - wh;
		    out[j*os] = f0 + wghp;
		    out[(j+n/4)*os] = f1 - I * wghm;
		    out[(j+n/2)*os] = f0 - wghp;
		    out[(j+3*n/4)*os] = f1 + I * wghm;
	       }

	       adds += (COMPLEX_MUL_ADDS * 2 + 2 * 6) * 2;
	       muls += (COMPLEX_MUL_MULS * 2) * 2;
	       twids += 2;
	  }

	  /* i == n/8 case (simpler multiply): */
	  if (i == n/8) {
	       Rt wabs = sqrt(0.5);
	       C f0 = out[i*os];
               C f1 = out[(i+n/4)*os];
               C wg = wabs * ((1.0 - I) * out[(i+n/2)*os]);
               C wh = wabs * ((1.0 + I) * out[(i+3*n/4)*os]);
               C wghp = wg + wh;
               C wghm = wg - wh;
               out[i*os] = f0 + wghp;
               out[(i+n/4)*os] = f1 - I * wghm;
               out[(i+n/2)*os] = f0 - wghp;
               out[(i+3*n/4)*os] = f1 + I * wghm;

	       adds += 2 * 2 + 2 * 6;
	       muls += 2 * 2;
	       twids += 1;
	  }
     }
}
