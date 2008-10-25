/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>
#include "kiss_fft.h"
#include "kiss_fftr.h"
#include "kiss_fftnd.h"
#include "kiss_fftndr.h"

BEGIN_BENCH_DOC
BENCH_DOC("name", "kissfft")
BENCH_DOC("version", "1.2.8")
BENCH_DOC("year", "2008")
BENCH_DOC("author", "Mark Borgerding")
BENCH_DOC("language", "C")
BENCH_DOC("url", "http://sourceforge.net/projects/kissfft/")
BENCH_DOC("copyright",
"Copyright (c) 2003-2006, Mark Borgerding\n"
"\n"
"All rights reserved.\n"
"\n"
"Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:\n"
"\n"
"    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.\n"
"    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.\n"
"    * Neither the author nor the names of any contributors may be used to endorse or promote products derived from this software without specific prior written permission.\n"
"\n"
	  "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return ((p->kind == PROBLEM_COMPLEX) ||
	     (p->rank > 0 && p->n[p->rank-1] % 2 == 0 
	      && p->kind == PROBLEM_REAL));
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_unpacked(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_unpacked(p, in, -1.0);
}

void copy_r2c(struct problem *p, bench_complex *out)
{
     copy_r2c_packed(p, out);
}

void copy_c2r(struct problem *p, bench_complex *in)
{
     copy_c2r_packed(p, in);
}

static void *cfg;

void setup(struct problem *p)
{
     size_t lenmem;
     if (p->kind == PROBLEM_COMPLEX) {
	  if (p->rank == 1) {
	       kiss_fft_alloc(p->n[0], p->sign == 1, 0, &lenmem);
	       cfg = bench_malloc(lenmem);
	       kiss_fft_alloc(p->n[0], p->sign == 1, cfg, &lenmem);
	  }
	  else {
	       unsigned int i;
	       int *dims;
	       dims = (int *) bench_malloc(sizeof(int) * p->rank);
	       for (i = 0; i < p->rank; ++i)
		    dims[i] = p->n[i];

	       kiss_fftnd_alloc(dims, p->rank, p->sign == 1, 0, &lenmem);
	       cfg = bench_malloc(lenmem);
	       kiss_fftnd_alloc(dims, p->rank, p->sign == 1, cfg, &lenmem);

	       bench_free(dims);
	  }
     }
     else if (p->rank == 1) {
	  kiss_fftr_alloc(p->n[0], p->sign == 1, 0, &lenmem);
	  cfg = bench_malloc(lenmem);
	  kiss_fftr_alloc(p->n[0], p->sign == 1, cfg, &lenmem);
     }
     else {
	       unsigned int i;
	       int *dims;
	       dims = (int *) bench_malloc(sizeof(int) * p->rank);
	       for (i = 0; i < p->rank; ++i)
		    dims[i] = p->n[i];

	       kiss_fftndr_alloc(dims, p->rank, p->sign == 1, 0, &lenmem);
	       cfg = bench_malloc(lenmem);
	       kiss_fftndr_alloc(dims, p->rank, p->sign == 1, cfg, &lenmem);

	       bench_free(dims);
     }
}

void doit(int iter, struct problem *p)
{
     int i;

     if (p->kind == PROBLEM_COMPLEX) {
	  kiss_fft_cpx *in = (kiss_fft_cpx *) p->in;
	  
	  if (p->in_place) {
	       if (p->rank == 1)
		    for (i = 0; i < iter; ++i) 
			 kiss_fft((kiss_fft_cfg) cfg, in, in);
	       else
		    for (i = 0; i < iter; ++i) 
			 kiss_fftnd((kiss_fftnd_cfg) cfg, in, in);
	  }
	  else {
	       kiss_fft_cpx *out = (kiss_fft_cpx *) p->out;
	       if (p->rank == 1)
		    for (i = 0; i < iter; ++i) 
			 kiss_fft((kiss_fft_cfg) cfg, in, out);
	       else
		    for (i = 0; i < iter; ++i) 
			 kiss_fftnd((kiss_fftnd_cfg) cfg, in, out);
	  }
     }
     else if (p->kind == PROBLEM_REAL) {
	  if (p->sign == 1) {
	       kiss_fft_cpx *in = (kiss_fft_cpx *) p->in;

	       if (p->in_place) {
		    if (p->rank == 1)
			 for (i = 0; i < iter; ++i) 
			      kiss_fftri((kiss_fftr_cfg) cfg, 
					 in, (kiss_fft_scalar *) in);
		    else
			 for (i = 0; i < iter; ++i) 
			      kiss_fftndri((kiss_fftndr_cfg) cfg, 
					   in, (kiss_fft_scalar *) in);
	       }
	       else {
		    kiss_fft_scalar *out = (kiss_fft_scalar *) p->out;
		    if (p->rank == 1)
			 for (i = 0; i < iter; ++i) 
			      kiss_fftri((kiss_fftr_cfg) cfg, in, out);
		    else
			 for (i = 0; i < iter; ++i) 
			      kiss_fftndri((kiss_fftndr_cfg) cfg, in, out);
	       }
	  }
	  else {
	       kiss_fft_scalar *in = (kiss_fft_scalar *) p->in;

	       if (p->in_place) {
		    if (p->rank == 1)
			 for (i = 0; i < iter; ++i) 
			      kiss_fftr((kiss_fftr_cfg) cfg, 
					in, (kiss_fft_cpx *) in);
		    else
			 for (i = 0; i < iter; ++i) 
			      kiss_fftndr((kiss_fftndr_cfg) cfg, 
					  in, (kiss_fft_cpx *) in);
	       }
	       else {
		    kiss_fft_cpx *out = (kiss_fft_cpx *) p->out;
		    if (p->rank == 1)
			 for (i = 0; i < iter; ++i) 
			      kiss_fftr((kiss_fftr_cfg) cfg, in, out);
		    else
			 for (i = 0; i < iter; ++i) 
			      kiss_fftndr((kiss_fftndr_cfg) cfg, in, out);
	       }
	  }
     }
}

void done(struct problem *p)
{
     bench_free(cfg);
     UNUSED(p);
}
