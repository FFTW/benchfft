/* this program is in the public domain */

#include "bench-user.h"
#include "complex.h"

BEGIN_BENCH_DOC
BENCH_DOC("name", "valkenburg")
BENCH_DOC("author", "Peter Valkenburg")
BENCH_DOC("year", "1987")
BENCH_DOC("language", "C")
BENCH_DOC("email", "p_valkenburg@hotmail.com")
BENCH_DOC("url", "http://www.jjj.de/fft/fft2.tgz")
BENCH_DOC("url-was-valid-on", "Mon Sep  2 18:10:40 EDT 2002")
BENCH_DOC("notes", "The forward transform has two extra passes to conjugate and scale input and output.")
BENCH_DOC("copyright",
"WARNING:  This package shows serious deficiencies if used in SDI-systems or "
"AEGIS-alikes.  So all of you defense-people: HANDS-OFF!")
END_BENCH_DOC

extern int fft (COMPLEX *in, unsigned int n, COMPLEX *out);
extern int rft (COMPLEX *in, unsigned int n, COMPLEX *out);
extern int W_init(unsigned int n);

int can_do(struct problem *p)
{
     return (DOUBLE_PRECISION &&
	     p->rank == 1 &&
	     p->kind == PROBLEM_COMPLEX &&
	     !problem_in_place(p));
}


void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
     W_init(p->n[0]);
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, -1);
}

void doit(int iter, struct problem *p)
{
     int n = p->n[0];
     COMPLEX *x = (COMPLEX *) p->in;
     COMPLEX *y = (COMPLEX *) p->out;
     int i;

     if (p->sign < 0)
	  for (i = 0; i < iter; ++i) {
	       fft(x, n, y);
	  }
     else
	  for (i = 0; i < iter; ++i) {
	       rft(x, n, y);
	  }
}

void done(struct problem *p)
{
     UNUSED(p);
}
