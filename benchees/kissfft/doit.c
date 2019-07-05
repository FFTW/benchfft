/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>
#include "kiss_fft.h"
#include "kiss_fftr.h"
#include "kiss_fftnd.h"
#include "kiss_fftndr.h"

BEGIN_BENCH_DOC
BENCH_DOC("name", "kissfft")
BENCH_DOC("version", "1.3.1")
BENCH_DOC("year", "2019")
BENCH_DOC("author", "Mark Borgerding")
BENCH_DOC("language", "C")
BENCH_DOC("url", "https://github.com/mborgerding/kissfft")
BENCH_DOC("copyright",
"Copyright (c) 2003-2010 Mark Borgerding . All rights reserved.\n"
"\n"
"KISS FFT is provided under:\n"
"\n"
"  SPDX-License-Identifier: BSD-3-Clause\n"
"\n"
"Being under the terms of the BSD 3-clause \"New\" or \"Revised\" License,\n"
"according with:\n"
"\n"
"  LICENSES/BSD-3-Clause\n"
)
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
