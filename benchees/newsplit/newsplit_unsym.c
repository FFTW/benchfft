/* Copyright (c) 2004 Massachusetts Institute of Technology
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

/* NOTE: This is essentially the same algorithm as newsplit.c 
   except with "unsymmetrized" scale factors.  These scale factors
   will have larger absolute magnitudes, and require larger twiddle
   tables, but should have the same number of flops and are easier
   to explain. */

/* scale factor pulled out of size n sub-DFT for mth output.
   By construction:
        scale(n, 0) = 1
	scale(n, m + n/4) = scale(n, m)
   and also, for 0 <= m <= n/8:
        omega_tan(n, m) = omega(n, m) * scale(n/4, m) / scale(n, m)
*/
Rt scale(int n, int m)
{
     if (n <= 4)
	  return 1.0;
     else {
	  int j = m % (n/4);
	  return cos((TWOPI / n) * j) * scale(n/4, j);
     }
}

void nsplit(int n, C *in0, C *in1, int is, C *out, int os)
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
	  nsplit(n/2, in0, in1 + is, is*2, out, os);
	  nsplitds(n/4, in1, in1 + 4*is, is*4, out + os * n/2, os);
	  nsplitds(n/4, in1 + (n-2)*is, in1 + 2*is, is*4, out + os * 3*n/4,os);
	  
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

	  for (i = 1; i < n/8; ++i) {
	       C w = omega(n, i) * scale(n/4, i);

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

	       adds += (COMPLEX_MUL_ADDS * 2 + 2 * 6);
	       muls += (COMPLEX_MUL_MULS * 2);
	  }

	  /* i == n/8 case (simpler multiply): */
	  if (i == n/8) {
	       Rt wabs = sqrt(0.5) * scale(n/4, i);
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
	  }

	  for (i = n/8+1; i < n/4; ++i) {
	       C w = omega(n, i) * scale(n/4, i);

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

	       adds += (COMPLEX_MUL_ADDS * 2 + 2 * 6);
	       muls += (COMPLEX_MUL_MULS * 2);
	  }

     }
}

/* as nsplit, but divide m-th output by scale(n, m) */
void nsplitds(int n, C *in0, C *in1, int is, C *out, int os)
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
	  nsplitds2(n/2, in0, in1 + is, is*2, out, os);
	  nsplitds(n/4, in1, in1 + 4*is, is*4, out + os * n/2, os);
	  nsplitds(n/4, in1 + (n-2)*is, in1 + 2*is, is*4, out + os * 3*n/4,os);
	  
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

	  for (i = 1; i < n/8; ++i) {
	       /* multiplying by omega_tan should only take
		  2 adds and 2 multiplies */
	       C w = omega_tan(n, i);
	       /*  = omega(n, i) * scale(n/4, i) / scale(n, i) */

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

	       adds += (2 * 2 + 2 * 6);
	       muls += (2 * 2);
	  }

	  /* i == n/8 case (simpler multiply): */
	  if (i == n/8) {
	       C f0 = out[i*os];
               C f1 = out[(i+n/4)*os];
               C wg = ((1.0 - I) * out[(i+n/2)*os]);
               C wh = ((1.0 + I) * out[(i+3*n/4)*os]);
               C wghp = wg + wh;
               C wghm = wg - wh;
               out[i*os] = f0 + wghp;
               out[(i+n/4)*os] = f1 - I * wghm;
               out[(i+n/2)*os] = f0 - wghp;
               out[(i+3*n/4)*os] = f1 + I * wghm;

	       adds += 2 * 2 + 2 * 6;
	  }

	  for (i = n/8+1; i < n/4; ++i) {
	       /* multiplying by omega_tan should only take
		  2 adds and 2 multiplies */
	       C w = omega_tan(n, i);
	       /*  = omega(n, i) * scale(n/4, i) / scale(n, i) */

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

	       adds += (2 * 2 + 2 * 6);
	       muls += (2 * 2);
	  }

     }
}

/* as nsplit, but divide m-th output by scale(2*n, m) */
/*                                       = wcos(2*n, m) * scale(n/2, m) */
void nsplitds2(int n, C *in0, C *in1, int is, C *out, int os)
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
	  nsplitds4(n/2, in0, in1 + is, is*2, out, os);
	  nsplitds(n/4, in1, in1 + 4*is, is*4, out + os * n/2, os);
	  nsplitds(n/4, in1 + (n-2)*is, in1 + 2*is, is*4, out + os * 3*n/4,os);
	  
	  /* i == 0 case (few multiplies): */
	  {
	       Rt s2 = 1.0 / scale(2*n, n/4);
	       C f0 = out[0*os];
               C f1 = out[(0+n/4)*os];
               C g = out[(0+n/2)*os];
               C h = out[(0+3*n/4)*os];
               C wghp = g + h;
               C wghm = s2 * (g - h);
               out[0*os] = f0 + wghp;
               out[(0+n/4)*os] = f1 - I * wghm;
               out[(0+n/2)*os] = f0 - wghp;
               out[(0+3*n/4)*os] = f1 + I * wghm;

	       adds += 2 * 6;
	       muls += 1 * 2;
	  }

	  for (i = 1; i < n/8; ++i) {
	       Ct wscale = wcos(n, i) * scale(n/4, i);
	       Ct w = omega_tan(n, i); /* *w should take 2 mults + 2 adds */

	       {
		    Rt s1 = 1.0 / scale(2*n, i);
		    Rt s2 = 1.0 / scale(2*n, i + n/4);
		    C f0 = out[i*os];
		    C f1 = out[(i+n/4)*os];
		    C wg = w * out[(i+n/2)*os];
		    C wh = conj(w) * out[(i+3*n/4)*os];
		    C wghp = (wg + wh) * (s1 * wscale);
		    C wghm = (wg - wh) * (s2 * wscale);
		    out[i*os] = f0 + wghp;
		    out[(i+n/4)*os] = f1 - I * wghm;
		    out[(i+n/2)*os] = f0 - wghp;
		    out[(i+3*n/4)*os] = f1 + I * wghm;
	       }
	       
	       adds += (2 * 2 + 2 * 6);
	       muls += (2 * 2 + 2 * 2);
	  }

	  /* i == n/8 case (simpler multiply): */
	  if (i == n/8) {
	       Rt s1 = 1.0 / scale(2*n, i);
	       Rt s2 = 1.0 / scale(2*n, i + n/4);
	       Rt wabs = sqrt(0.5) * scale(n/4, i);
	       C f0 = out[i*os];
               C f1 = out[(i+n/4)*os];
               C wg = ((1.0 - I) * out[(i+n/2)*os]);
               C wh = ((1.0 + I) * out[(i+3*n/4)*os]);
               C wghp = (wg + wh) * (wabs * s1);
               C wghm = (wg - wh) * (wabs * s2);
               out[i*os] = f0 + wghp;
               out[(i+n/4)*os] = f1 - I * wghm;
               out[(i+n/2)*os] = f0 - wghp;
               out[(i+3*n/4)*os] = f1 + I * wghm;

	       adds += 2 * 2 + 2 * 6;
	       muls += 2 * 2;
	  }

	  for (i = n/8+1; i < n/4; ++i) {
	       Ct wscale = wcos(n, i) * scale(n/4, i);
	       Ct w = omega_tan(n, i); /* *w should take 2 mults + 2 adds */

	       {
		    Rt s1 = 1.0 / scale(2*n, i);
		    Rt s2 = 1.0 / scale(2*n, i + n/4);
		    C f0 = out[i*os];
		    C f1 = out[(i+n/4)*os];
		    C wg = w * out[(i+n/2)*os];
		    C wh = conj(w) * out[(i+3*n/4)*os];
		    C wghp = (wg + wh) * (s1 * wscale);
		    C wghm = (wg - wh) * (s2 * wscale);
		    out[i*os] = f0 + wghp;
		    out[(i+n/4)*os] = f1 - I * wghm;
		    out[(i+n/2)*os] = f0 - wghp;
		    out[(i+3*n/4)*os] = f1 + I * wghm;
	       }
	       
	       adds += (2 * 2 + 2 * 6);
	       muls += (2 * 2 + 2 * 2);
	  }

     }
}

/* as nsplit, but divide m-th output by scale(4*n, m) */
void nsplitds4(int n, C *in0, C *in1, int is, C *out, int os)
{
     if (n == 1)
	  out[0] = *in0;
     else if (n == 2) {
	  Rt s1 = sqrt(2.0); /* = 1.0 / scale(4*n, 1); */
	  C a, b;
	  a = *in0;
	  b = *in1;
	  out[0] = a + b;
	  out[os] = (a - b) * s1;

	  adds += 2 * 2;
	  muls += 1 * 2;
     }
     else {
	  int i;
	  nsplitds2(n/2, in0, in1 + is, is*2, out, os);
	  nsplitds(n/4, in1, in1 + 4*is, is*4, out + os * n/2, os);
	  nsplitds(n/4, in1 + (n-2)*is, in1 + 2*is, is*4, out + os * 3*n/4,os);
	  
	  /* i == 0 case (few multiplies): */
	  {
	       Rt s2 = 1.0 / scale(4*n, n/4);
	       Rt s3 = 1.0 / scale(4*n, n/2);
	       Rt s4 = 1.0 / scale(4*n, 3*n/4);
	       C f0 = out[0*os];
               C f1 = out[(0+n/4)*os];
               C g = out[(0+n/2)*os];
               C h = out[(0+3*n/4)*os];
               C wghp = g + h;
               C wghm = (g - h);
               out[0*os] = f0 + wghp;
               out[(0+n/4)*os] = (f1 - I * wghm) * s2;
               out[(0+n/2)*os] = (f0 - wghp) * s3;
               out[(0+3*n/4)*os] = (f1 + I * wghm) * s4;

	       adds += 2 * 6;
	       muls += 2 * 3;
	  }

	  for (i = 1; i < n/8; ++i) {
	       Ct wscale = wcos(n, i) * scale(n/4, i);
	       Ct w = omega_tan(n, i); /* *w should take 2 mults + 2 adds */

	       {
		    Rt s1 = 1.0 / scale(4*n, i);
		    Rt s2 = 1.0 / scale(4*n, i + n/4);
		    Rt s3 = 1.0 / scale(4*n, i + n/2);
		    Rt s4 = 1.0 / scale(4*n, i + 3*n/4);
		    C f0 = out[i*os];
		    C f1 = out[(i+n/4)*os];
		    C wg = w * out[(i+n/2)*os];
		    C wh = conj(w) * out[(i+3*n/4)*os];
		    C wghp = wg + wh;
		    C wghm = wg - wh;
		    out[i*os] = (f0 + wghp) * (s1 * wscale);
		    out[(i+n/4)*os] = (f1 - I * wghm) * (s2 * wscale);
		    out[(i+n/2)*os] = (f0 - wghp) * (s3 * wscale);
		    out[(i+3*n/4)*os] = (f1 + I * wghm) * (s4 * wscale);
	       }
	       
	       adds += (2 * 2 + 2 * 6);
	       muls += (2 * 2 + 2 * 4);
	  }

	  /* i == n/8 case (simpler multiply): */
	  if (i == n/8) {
	       Rt s1 = 1.0 / scale(4*n, i);
	       Rt s2 = 1.0 / scale(4*n, i + n/4);
	       Rt s3 = 1.0 / scale(4*n, i + n/2);
	       Rt s4 = 1.0 / scale(4*n, i + 3*n/4);
	       Rt wabs = sqrt(0.5) * scale(n/4, i);
	       C f0 = out[i*os];
               C f1 = out[(i+n/4)*os];
               C wg = ((1.0 - I) * out[(i+n/2)*os]);
               C wh = ((1.0 + I) * out[(i+3*n/4)*os]);
               C wghp = wg + wh;
               C wghm = wg - wh;
               out[i*os] = (f0 + wghp) * (s1 * wabs);
               out[(i+n/4)*os] = (f1 - I * wghm) * (s2 * wabs);
               out[(i+n/2)*os] = (f0 - wghp) * (s3 * wabs);
               out[(i+3*n/4)*os] = (f1 + I * wghm) * (s4 * wabs);

	       adds += 2 * 2 + 2 * 6;
	       muls += 2 * 4;
	  }

	  for (i = n/8+1; i < n/4; ++i) {
	       Ct wscale = wcos(n, i) * scale(n/4, i);
	       Ct w = omega_tan(n, i); /* *w should take 2 mults + 2 adds */

	       {
		    Rt s1 = 1.0 / scale(4*n, i);
		    Rt s2 = 1.0 / scale(4*n, i + n/4);
		    Rt s3 = 1.0 / scale(4*n, i + n/2);
		    Rt s4 = 1.0 / scale(4*n, i + 3*n/4);
		    C f0 = out[i*os];
		    C f1 = out[(i+n/4)*os];
		    C wg = w * out[(i+n/2)*os];
		    C wh = conj(w) * out[(i+3*n/4)*os];
		    C wghp = wg + wh;
		    C wghm = wg - wh;
		    out[i*os] = (f0 + wghp) * (s1 * wscale);
		    out[(i+n/4)*os] = (f1 - I * wghm) * (s2 * wscale);
		    out[(i+n/2)*os] = (f0 - wghp) * (s3 * wscale);
		    out[(i+3*n/4)*os] = (f1 + I * wghm) * (s4 * wscale);
	       }
	       
	       adds += (2 * 2 + 2 * 6);
	       muls += (2 * 2 + 2 * 4);
	  }

     }
}
