/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "emayer")
BENCH_DOC("author", "Ernst W. Mayer")
BENCH_DOC("year", "1997")
BENCH_DOC("language", "Fortran 90")
BENCH_DOC("email", "ewmayer@aol.com")
BENCH_DOC("notes", "A later (radix-16) version of this code was published in the Mersenne Prime Freeware collection; the later version merges forward and backwards transforms and so is more difficult to benchmark.")
BENCH_DOC("notes", "The original URL for this code is now defunct.")
BENCH_DOC("copyright",
"Copyright 1997 by Ernst W. Mayer. This program may be used and\n"
"redistributed freely as long as this header is included. You may\n"
"modify this program freely, as long as any redistribution contains the\n"
"original header and a summary of any changes made.\n"
"\n"
"This software is offered *as is*, without warranty of any kind. It is\n"
"intended only for the private recreational use of our audience. If it\n"
"causes your CPU to go *poof*, that's tough.\n"
"\n"
"Any modifications of this program intended for distribution MUST be\n"
"put into the public domain, or provided freely to the original\n"
"author. Under no circumstances is this program, or any program derived\n"
"from it or making use of it in any way, to be used for commercial\n"
"activities of any kind without the express written consent of the\n"
"original author. (Read: you aren't allowed to make money off it\n"
"without my say-so.)\n"
"\n"
"If you use this program or any program derived from it in work leading\n"
"to a publication, proper acknowledgement should be made in said\n"
"publication.\n"
"\n"
"The author appreciates if users send bug reports, inquiries related to\n"
"the software, or reports of useful modifications or novel and\n"
"interesting applications to him at the following address:\n"
"\n"
"Ernst Mayer\n"
"Dept. of Mechanical and Aerospace Engineering\n"
"Case Western Reserve University\n"
"10900 Euclid Avenue\n"
"Cleveland, OH 44106-7222 USA\n"
"E-mail: mayer@nigel.mae.cwru.edu or ewmayer@aol.com\n")
END_BENCH_DOC


int can_do(struct problem *p)
{
     return (DOUBLE_PRECISION &&
	     p->rank == 1 &&
	     problem_power_of_two(p, 1) &&
	     ((p->kind == PROBLEM_COMPLEX && p->n[0] >=8 && p->n[0] <=(1<<18))
	      ||
	      (p->kind == PROBLEM_REAL && p->n[0] >= 16 && p->n[0] <= (1<<19)))
	  );
}

#define FFT F77_FUNC_(dlanczos_fwd,DLANCZOS_FWD)
extern void FFT(bench_complex *a, int *two_n);

#define IFFT F77_FUNC_(dlanczos_rev,DLANCZOS_REV)
extern void IFFT(bench_complex *a, int *two_n);

#define REAL_FFT F77_FUNC_(dfft_fwd,DFFT_FWD)
extern void REAL_FFT(bench_real *a, int *n);

#define REAL_IFFT F77_FUNC_(dfft_rev,DFFT_REV)
extern void REAL_IFFT(bench_real *a, int *n);

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));

     /* Call FFT once to initialize things before benchmarking: */
     doit(1, p);
}

void doit(int iter, struct problem *p)
{
     int i;
     if (p->kind == PROBLEM_COMPLEX) {
	  bench_complex *in = (bench_complex *) p->in;
	  int two_n = 2 * p->n[0];

	  if (p->sign > 0)
	       for (i = 0; i < iter; ++i)
		    FFT(in, &two_n);
	  else
	       for (i = 0; i < iter; ++i)
		    IFFT(in, &two_n);
     }
     else /* PROBLEM_REAL */ {
	  bench_real *in = (bench_real *) p->in;
	  int n = p->n[0];

	  if (p->sign < 0)
	       for (i = 0; i < iter; ++i)
		    REAL_FFT(in, &n);
	  else
	       for (i = 0; i < iter; ++i)
		    REAL_IFFT(in, &n);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}

/* Backwards real transform is scaled by 0.5 */
void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     if (p->kind == PROBLEM_REAL && p->sign > 0) {
	  bench_complex x;
	  c_re(x) = 2.0;
	  c_im(x) = 0;
	  
	  cascale(out, p->size, x);
     }
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_1d_packed(p, out, 1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_packed(p, in, 1.0);
}

