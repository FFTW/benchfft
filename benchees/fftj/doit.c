/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
#if defined(FFTJ_FORTRAN)
BENCH_DOC("name", "fftj-fortran")
BENCH_DOC("language", "Fortran 77")
#elif defined (FFTJ_SSE)
BENCH_DOC("name", "fftj-sse")
BENCH_DOC("language", "Fortran 77, x86 assembly")
#else
#error "unknown fftj version"
#endif
BENCH_DOC("author", "Keiichi Ishioka")
BENCH_DOC("email", " ishioka@gfd-dennou.org")
BENCH_DOC("year", "2008")
END_BENCH_DOC

#define MAXN 1024

#define FJCINI F77_FUNC(fjcini,FJCINI)
#define FJCRUN F77_FUNC(fjcrun,FJCRUN)

extern void FJCINI(int *N, int *IO, int *IBF, void *IP[1],
		   bench_complex *ZT);

extern void FJCRUN(bench_complex *ZIN,
		   bench_complex *ZOUT,
		   bench_complex *ZWORK,
		   bench_complex *ZT,
		   void *IP[1]);

int can_do(struct problem *p)
{
     return (p->rank == 1 && 
	     p->kind == PROBLEM_COMPLEX &&
	     power_of_two(p->n[0]) &&
	     p->n[0] <= MAXN);
}

bench_complex *ZT, *ZWORK;
void *IP[1];

void setup(struct problem *p)
{
     int N = p->n[0];
     int IO = p->in == p->out ? 1 : 2;
     int IBF = p->sign == -1 ? 2 : 1;

     ZT = bench_malloc(2 * MAXN * sizeof(*ZT));
     ZWORK = bench_malloc(MAXN * sizeof(*ZWORK));
     FJCINI(&N, &IO, &IBF, IP, ZT);
}

void doit(int iter, struct problem *p)
{
     int i;
     bench_complex *in = p->in;
     bench_complex *out = p->out;

     if (in == out) {
	  for (i = 0; i < iter; ++i) {
	       FJCRUN(in, ZWORK, ZWORK, ZT, IP);
	  }
     } else {
	  for (i = 0; i < iter; ++i) {
	       FJCRUN(in, out, ZWORK, ZT, IP);
	  }
     }
}

void done(struct problem *p)
{
     bench_free(ZWORK);
     bench_free(ZT);
     UNUSED(p);
}
