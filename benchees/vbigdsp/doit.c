/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

#include "vBigDSP.h"

BEGIN_BENCH_DOC
BENCH_DOC("name", "vbigdsp")
BENCH_DOC("author", "Apple Computer, Inc.")
BENCH_DOC("author", "Greg Allen")
BENCH_DOC("email", "gallen@arlut.utexas.edu")
BENCH_DOC("version", "20030322")
BENCH_DOC("year", "2003")
BENCH_DOC("copyright",
"Copyright (c) 1999 by Apple Computer, Inc.  All rights reserved.\n"
"\n"
"IMPORTANT: This Apple software is supplied to you by Apple Computer,\n"
"Inc.  ('Apple') in consideration of your agreement to the following\n"
"terms, and your use, installation, modification or redistribution of\n"
"this Apple software constitutes acceptance of these terms.  If you do\n"
"not agree with these terms, please do not use, install, modify or\n"
"redistribute this Apple software.\n"
"\n"
"In consideration of your agreement to abide by the following terms,\n"
"and subject to these terms, Apple grants you a personal, non-exclusive\n"
"license, under Apple's copyrights in this original Apple software (the\n"
"'Apple Software'), to use, reproduce, modify and redistribute the\n"
"Apple Software, with or without modifications, in source and/or binary\n"
"forms; provided that if you redistribute the Apple Software in its\n"
"entirety and without modifications, you must retain this notice and\n"
"the following text and disclaimers in all such redistributions of the\n"
"Apple Software.  Neither the name, trademarks, service marks or logos\n"
"of Apple Computer, Inc. may be used to endorse or promote products\n"
"derived from the Apple Software without specific prior written\n"
"permission from Apple.  Except as expressly stated in this notice, no\n"
"other rights or licenses, express or implied, are granted by Apple\n"
"herein, including but not limited to any patent rights that may be\n"
"infringed by your derivative works or by other works in which the\n"
"Apple Software may be incorporated.\n"
"\n"
"The Apple Software is provided by Apple on an 'AS IS' basis.  APPLE\n"
"MAKES NO WARRANTIES, EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION\n"
"THE IMPLIED WARRANTIES OF NON-INFRINGEMENT, MERCHANTABILITY AND\n"
"FITNESS FOR A PARTICULAR PURPOSE, REGARDING THE APPLE SOFTWARE OR ITS\n"
"USE AND OPERATION ALONE OR IN COMBINATION WITH YOUR PRODUCTS.  IN NO\n"
"EVENT SHALL APPLE BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL OR\n"
"CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF\n"
"SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR\n"
"BUSINESS INTERRUPTION) ARISING IN ANY WAY OUT OF THE USE,\n"
"REPRODUCTION, MODIFICATION AND/OR DISTRIBUTION OF THE APPLE SOFTWARE,\n"
"HOWEVER CAUSED AND WHETHER UNDER THEORY OF CONTRACT, TORT (INCLUDING\n"
"NEGLIGENCE), STRICT LIABILITY OR OTHERWISE, EVEN IF APPLE HAS BEEN\n"
"ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n")
BENCH_DOC("notes", "Sample code provided by Apple, modified by Greg Allen to work on Linux/PowerPC.")
BENCH_DOC("notes",
	  "Apple's description: Though arbitrary signal lengths (i.e. all powers of 2) are handled, our design emphasis is on very long signals (length N >= 2^16 and on into the millions), for which cache considerations are paramount. The core of the library is a particular variant of full-complex FFT that for signal length N = 2^10 executes at 1.15 giga ops (500 MHz G4). This cache-friendly, core FFT plays a dominant role in the long-signal cases such as two-dimensional FFT and convolution. More important perhaps than the core performance benchmark is the manner in which one can sift through the myriad prevailing (and new) FFT frameworks, to arrive at a suitable such framework for the Velocity Engine.")
BENCH_DOC("url", "http://findsabrina.org/altivec/")
BENCH_DOC("url", "http://developer.apple.com/samplecode/Sample_Code/Devices_and_Hardware/Velocity_Engine/VelEng_FFT.htm")
BENCH_DOC("url-was-valid-on", "Thu Mar 27 20:18:01 EST 2003")
BENCH_DOC("bibitem",
	  "R. Crandall and J. Klivington, Supercomputer-style FFT library "
	  "for Apple G4, Apple Technical Report (Jan. 2000).")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (SINGLE_PRECISION &&
	     p->rank <= 2 &&
	     p->n[0] >= 2 &&
	     (p->kind == PROBLEM_COMPLEX 
	      || p->rank == 1
	      || (p->n[0] >= 8 && p->n[1] >= 8)) &&
             problem_power_of_two(p, 1));
}

static int imax2(int a, int b) { return (a > b ? a : b); }
static int imin2(int a, int b) { return (a < b ? a : b); }

unsigned long width = 0, height = 0;
long length = 0;
long isign = 0;

void setup(struct problem *p)
{
     isign = p->sign < 0 ? IFLAG_FFT_FORWARD : IFLAG_FFT_INVERSE;
     length = p->size;
     if (p->rank > 1) {
	  width = p->n[1];
	  height = p->n[0];
     }
     doit(1, p); /* fft tables are initialized on first call */
}

void doit(int iter, struct problem *p)
{
     float *data = (float *) p->in;
     int i;

     if (p->kind == PROBLEM_COMPLEX) {
	  switch (p->rank) {
	      case 1:
		   for (i = 0; i < iter; ++i) {
			FFTComplex(data, length, isign);
		   }
		   break;
	      case 2:
		   for (i = 0; i < iter; ++i) {
			FFT2DComplex(data, width, height, isign);
		   }
		   break;
	  }
     }
     else { /* PROBLEM_REAL */
	  if (isign == IFLAG_FFT_FORWARD) {
	       switch (p->rank) {
		   case 1: 
			for (i = 0; i < iter; ++i) {
			     FFTRealForward(data, length);
			}
			break;
		   case 2: 
			for (i = 0; i < iter; ++i) {
			     FFT2DRealForward(data, width, height);
			}
			break;
	       }
	  }
	  else { /* BACKWARD */
	       switch (p->rank) {
		   case 1: 
			for (i = 0; i < iter; ++i) {
			     FFTRealInverse(data, length);
			}
			break;
		   case 2: 
			for (i = 0; i < iter; ++i) {
			     FFT2DRealInverse(data, width, height);
			}
			break;
	       }
	  }
     }
}

void done(struct problem *p)
{
     ShutdownFFT();
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, 1);
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     if (p->rank == 1) {
	  copy_h2c_1d_packed(p, out, -1);
     }
     else { /* insane format */
	  bench_real *data = (float *) p->in;
	  int H = height, W = width;
	  int i, j;
	  
	  c_re(out[0]) = data[0];
	  c_im(out[0]) = 0;

	  c_re(out[H/2 * W]) = data[W];
	  c_im(out[H/2 * W]) = 0;

	  /* rest of first column */
	  for (i = 1; i < H/2; ++i) {
	       c_re(out[(H-i)*W]) = c_re(out[i*W]) = data[2*i*W];
	       c_im(out[(H-i)*W]) = -(c_im(out[i*W]) = data[(2*i+1)*W]);
	  }

	  c_re(out[W/2]) = data[1];
	  c_im(out[W/2]) = 0;

	  c_re(out[H/2 * W + W/2]) = data[W + 1];
	  c_im(out[H/2 * W + W/2]) = 0;

	  /* rest of middle column */
	  for (i = 1; i < H/2; ++i) {
	       c_re(out[(H-i)*W+W/2]) = c_re(out[i*W+W/2]) = data[2*i*W + 1];
	       c_im(out[(H-i)*W+W/2]) = -(c_im(out[i*W+W/2]) 
					  = data[(2*i+1)*W + 1]);
	  }

	  /* rest of columns */
	  for (i = 1; i < H; ++i) {
	       for (j = 1; j < W/2; ++j) {
		    c_re(out[(H-i)*W + (W-j)]) 
			 = c_re(out[i*W + j]) = data[i*W + 2*j];
		    c_im(out[(H-i)*W + (W-j)]) 
			 = -(c_im(out[i*W + j]) = data[i*W + 2*j+1]);
	       }
	  }

	  /* rest of first row */
	  for (j = 1; j < W/2; ++j) {
	       c_re(out[(W-j)]) = c_re(out[j]) = data[2*j];
	       c_im(out[(W-j)]) = -(c_im(out[j]) = data[2*j+1]);
	  }
     }
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     if (p->rank == 1) {
	  copy_c2h_1d_packed(p, in, -1);
     }
     else { /* insane format */
	  bench_real *data = (float *) p->in;
	  int H = height, W = width;
	  int i, j;
	  
	  data[0] = c_re(in[0]);

	  data[W] = c_re(in[H/2 * W]);

	  /* rest of first column */
	  for (i = 1; i < H/2; ++i) {
	       data[2*i*W] = c_re(in[i*W]);
	       data[(2*i+1)*W] = c_im(in[i*W]);
	  }
	  
	  data[1] = c_re(in[W/2]);

	  data[W + 1] = c_re(in[H/2 * W + W/2]);

	  /* rest of middle column */
	  for (i = 1; i < H/2; ++i) {
	       data[2*i*W + 1] = c_re(in[i*W+W/2]);
	       data[(2*i+1)*W + 1] = c_im(in[i*W+W/2]);
	  }

	  /* rest of columns */
	  for (i = 1; i < H; ++i) {
	       for (j = 1; j < W/2; ++j) {
		    data[i*W + 2*j] = c_re(in[i*W + j]);
		    data[i*W + 2*j+1] = c_im(in[i*W + j]);
	       }
	  }

	  /* rest of first row */
	  for (j = 1; j < W/2; ++j) {
	       data[2*j] = c_re(in[j]);
	       data[2*j+1] = c_im(in[j]);
	  }
     }
}
