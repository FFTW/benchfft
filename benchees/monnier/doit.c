/* this program is in the public domain */

#include "bench-user.h"

BEGIN_BENCH_DOC
BENCH_DOC("name", "monnier")
BENCH_DOC("author", "Yves Monnier")
BENCH_DOC("year", "1995")
BENCH_DOC("language", "C")
BENCH_DOC("url", "http://www.igd.u-bordeaux.fr/~monnier/fft.html")
BENCH_DOC("url-was-valid-on", "August 1998 (?)")
BENCH_DOC("notes", "Slightly modified for original benchfft; "
	  "original version has disappeared from the Web.")
BENCH_DOC("notes", "Although a pre-computed trigonometric table is used, "
	  "the phase angle in the table is calculated by repeated addition, "
	  "resulting in poor accuracy even for moderate transform sizes.")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (p->rank == 1 &&
	     p->kind == PROBLEM_COMPLEX &&
	     problem_in_place(p));
}
	  
int ouvre_fft(int taille, bench_real sens);
int ferme_fft(void);
int fft(bench_real *donnees);

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, -1);
}

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
     BENCH_ASSERT(!ouvre_fft(p->n[0], p->sign));
}

void doit(int iter, struct problem *p)
{
     bench_real *x = (bench_real *) p->in;
     int i;

     for (i = 0; i < iter; ++i)
	  fft(x);
}

void done(struct problem *p)
{
     UNUSED(p);
     BENCH_ASSERT(!ferme_fft());
}
