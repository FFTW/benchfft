/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "green-ffts-2.0")
BENCH_DOC("author", "John Green")
BENCH_DOC("email", "green_jt@vsdec.npt.nuwc.navy.mil")
BENCH_DOC("language", "C")
BENCH_DOC("notes", "The backward transform is scaled")
END_BENCH_DOC

#include "fftext.h"

static int M;

int can_do(struct problem *p)
{
     return (sizeof(bench_real) == sizeof(float) &&
	     problem_complex_power_of_two(p, 1));
}

void problem_ccopy_to(struct problem *p, bench_complex *out)
{
     cacopy(p->out, out, p->size);
     if (p->sign == 1) { 	  /* undo the scaling */
	  bench_complex x;
	  c_re(x) = p->size;
	  c_im(x) = 0;

	  cascale(out, p->size, x);
     }
}

void setup(struct problem *p)
{
     int n = p->n[0];

     M = 0;

     while (n > 1) {
	  M += 1;
	  n /= 2;
     }

     fftInit(M);
}

void doit(int iter, struct problem *p)
{
     int i;
     int n = p->n[0];
     void *in = p->in;
     int m = M; /* cache */

     if (p->sign == -1) {
	  for (i = 0; i < iter; ++i) {
	       ffts(in, m, 1);
	  }
     } else {
	  for (i = 0; i < iter; ++i) {
	       iffts(in, m, 1);
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     fftFree();
}
