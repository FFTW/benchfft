/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "morris82")
BENCH_DOC("author", "L. Robert Morris")
BENCH_DOC("year", "1982")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("notes", "Received in personal communication from the author, Apr 21 2003")
BENCH_DOC("notes", "Complex data are stored in separate real/imag arrays.")
BENCH_DOC("notes", "Published in: Digital Signal Processing Software, DSPS Inc, 1983")
END_BENCH_DOC

extern void F77_FUNC(fft4, FFT4)();
extern void F77_FUNC(inir4, INIR4)();
static int m;

void copy_c2c_from(struct problem *p, bench_complex *in)
{
     bench_real *x = (bench_real *) p->in;
     /* 1024 matches magic number in morris.f */
     copy_c2ri(in, x, x + 1024, p->size);
}

void copy_c2c_to(struct problem *p, bench_complex *out)
{
     bench_real *x = (bench_real *) p->out;
     /* 1024 matches magic number in morris.f */
     copy_ri2c(x, x + 1024, out, p->size);
}

void problem_alloc(struct problem *p)
{
     /* the routine expects arguments in common block AA */
     extern char F77_FUNC(aa, AA)[];

     p->in = p->out = F77_FUNC(aa, AA);
}


void problem_free(struct problem *p)
{
     UNUSED(p);
     /* nothing to deallocate */
}

int can_do(struct problem *p)
{
     return (SINGLE_PRECISION &&
	     p->rank == 1 &&
	     p->sign == -1 &&
	     problem_complex_power_of_two(p, 1) &&
	     ((log_2(p->n[0]) & 1) == 0) && 	     /* power of 4 */
	     p->n[0] > 4 &&
	     p->n[0] <= 1024);
}


void setup(struct problem *p)
{
     m = log_2(p->n[0]) / 2;
     BENCH_ASSERT(can_do(p));
     F77_FUNC(inir4, INIR4)(&m);
}

void doit(int iter, struct problem *p)
{
     int i;

     for (i = 0; i < iter; ++i) {
	  F77_FUNC(fft4, FFT4)(&p->n[0], &m);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
