/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "dsp79-singleton")
BENCH_DOC("package", "Programs for Digital Signal Processing")
BENCH_DOC("author", "Richard C. Singleton")
BENCH_DOC("year", "1979, derived from Singleton 1968")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("notes", "DSP Section 1.4")
BENCH_DOC("bibitem", 
	  "R. C. Singleton, An algorithm for computing the mixed radix fast Fourier transform, IEEE Trans. on Audio and Electroacoustics AU-17, no. 2, p. 93-103 (June, 1969).")
BENCH_DOC("notes", "The backward transform is scaled")
BENCH_DOC("bibitem", "IEEE DSP Committee, Programs for Digital Signal Processing (IEEE Press, 1979). ISBN 0-471-05962-5. (Out of Print)")
END_BENCH_DOC


void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, 1);
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_unpacked(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_unpacked(p, in, -1.0);
}

int can_do(struct problem *p)
{
     return (SINGLE_PRECISION &&
	     p->rank == 1 &&
	     /* real transform works only for even N */
	     (p->kind == PROBLEM_COMPLEX || (p->n[0] & 1) == 0) &&
	     problem_in_place(p) &&
p->n[0] != 1960 &&   /* doesn't like this */
	     /* not quite... also works for some N larger than 12754584 */
	     p->n[0] < 12754584); 
}

extern void F77_FUNC(fft, FFT)();
extern void F77_FUNC(reals, REALS)();

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
}

void doit(int iter, struct problem *p)
{
     int i;
     int nseg = 1;
     int nspn = 1;

     if (p->kind == PROBLEM_COMPLEX) {
	  int n = p->n[0];
	  int isn = (p->sign > 0) ? 2 : -2;
	  bench_complex *in = p->in;

	  for (i = 0; i < iter; ++i) {
	       F77_FUNC(fft, FFT)(&c_re(in[0]), &c_im(in[0]),
				  &nseg, &n, &nspn, &isn);
	  }
     } else {
	  int n = p->n[0] / 2;
	  bench_real *in = p->in;

	  if (p->sign == -1) {
	       int isn = -2;
	       for (i = 0; i < iter; ++i) {
		    /* calling sequence suggested by source code: */
		    F77_FUNC(fft, FFT)(in, in + 1, &nseg, &n, &nspn, &isn);
		    F77_FUNC(reals, REALS)(in, in + 1, &n, &isn);
	       }
	  } else {
	       int isn = 2;
	       for (i = 0; i < iter; ++i) {
		    /* calling sequence suggested by source code: */
		    F77_FUNC(reals, REALS)(in, in + 1, &n, &isn);
		    F77_FUNC(fft, FFT)(in, in + 1, &nseg, &n, &nspn, &isn);
	       }
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
