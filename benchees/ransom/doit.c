/* this program is in the public domain */

#include "bench-user.h"

BEGIN_BENCH_DOC
BENCH_DOC("name", "ransom")
BENCH_DOC("author", "Scott Ransom")
BENCH_DOC("email", "ransom@physics.mcgill.ca")
BENCH_DOC("year", "1997")
BENCH_DOC("language", "C")
BENCH_DOC("notes", "Received in personal communication with the author.")
BENCH_DOC("notes", "Uses a six-step FFT algorithm.  Derived from Mikko Tommila's Apfloat arbitrary-precision arithmetic package, which performed a six-step number-theoretic transform.")
BENCH_DOC("notes", "Backwards transform works by two extra passes to conjugate input and output.")
BENCH_DOC("copyright",
"         Copyright (c) 1997 Scott M. Ransom\n"
"\n"
"         *** LEGAL NOTICE: ***\n"
"         (borrowed and modified from Joerg Arndt's fxt package)\n"
"\n"
"1.) The code herein is freeware.  You may use it at no cost\n"
"    and give it away to other people.\n"
"\n"
"2.) You should always attach a copy of the original package\n"
"    to your derived work that you give away (at least a pointer to it).\n"
"\n"
"3.) You are NOT allowed to make money with this code in any way.\n"
"\n"
"4.) I make no guarantee about the code.  It could, and probably does,\n"
"    contain some nasty bug somewhere.\n"
"\n"
"5.) Be nice and let me know if you use this package for something\n"
"    interesting/noticable/scientific.\n"
"          *** End of Legal Notice ***\n")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (SINGLE_PRECISION &&
	     p->rank == 1 &&
	     problem_power_of_two(p, 1) &&
	     p->kind == PROBLEM_COMPLEX && /* real transform is buggy */
	     p->n[0] > (p->kind == PROBLEM_COMPLEX ? 2U : 4U));
}

extern void tablesixstepfft(float *indata, long nn, int isign);
extern void realfft(float idata[], unsigned long n, int isign);

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_1d_packed(p, out, -1.0);

     /* correct for an apparent bug in the package */
     c_im(out[p->n[0]/4]) *= -1;
     c_im(out[p->n[0] - p->n[0]/4]) *= -1;
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_packed(p, in, -1.0);

     { /* correct for an apparent bug in the package */
	  bench_complex *x = (bench_complex *) p->in;
	  c_im(x[p->n[0]/4]) *= -1;
     }
     /* this doesn't seem to fix all the bugs in the inverse real transform,
	though.  punt. */
}

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     if (p->kind == PROBLEM_REAL && p->sign > 0) {
          bench_complex x;
          c_re(x) = 2;
          c_im(x) = 0;
          cascale(out, p->size, x);
     }
}

void doit(int iter, struct problem *p)
{
     int n = p->n[0], isign = p->sign;
     bench_real *x = (bench_real *) p->in;
     int i;

     if (p->kind == PROBLEM_COMPLEX)
	  for (i = 0; i < iter; ++i) {
	       tablesixstepfft(x, n, isign);
	  }
     else {
	  for (i = 0; i < iter; ++i) {
	       realfft(x, n, isign);
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
