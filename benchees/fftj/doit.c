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
BENCH_DOC("url", "http://www.gfd-dennou.org/arch/ishioka/fftj/")
BENCH_DOC("url-was-valid-on", "Fri Oct 24 11:09:09 EDT 2008")
END_BENCH_DOC

#define FJCINI F77_FUNC(fjcini,FJCINI)
#define FJCRUN F77_FUNC(fjcrun,FJCRUN)
#define FJRINI F77_FUNC(fjrini,FJRINI)
#define FJRRUN F77_FUNC(fjrrun,FJRRUN)

extern void FJCINI(int *N, int *IO, int *IBF, void *IP[1],
		   bench_complex *ZT);

extern void FJCRUN(bench_complex *ZIN,
		   bench_complex *ZOUT,
		   bench_complex *ZWORK,
		   bench_complex *ZT,
		   void *IP[1]);

extern void FJRINI(int *N, int *IO, int *IBF, void *IP[1],
		   bench_real *XT);

extern void FJRRUN(bench_real *XIN,
		   bench_real *XOUT,
		   bench_real *XWORK,
		   bench_real *XT,
		   void *IP[1]);

int can_do(struct problem *p)
{
     return (p->rank == 1 && 
	     power_of_two(p->n[0]) &&
	     p->n[0] > 1 &&
	     ((p->kind == PROBLEM_COMPLEX && p->n[0] <= 1024) ||
	      (p->kind == PROBLEM_REAL && p->n[0] <= 2048)));
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_1d_packed(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_packed(p, in, -1.0);
}

void copy_r2c(struct problem *p, bench_complex *out)
{
     copy_r2c_packed(p, out);
}

void copy_c2r(struct problem *p, bench_complex *in)
{
     copy_c2r_packed(p, in);
}


bench_complex *ZT, *ZWORK;
bench_real *XT, *XWORK;
void *IP[1];

void setup(struct problem *p)
{
     int N = p->n[0];
     int IO = p->in == p->out ? 1 : 2;
     int IBF = p->sign == -1 ? 2 : 1;

     if (p->kind == PROBLEM_COMPLEX) {
	  ZT = bench_malloc(2 * N * sizeof(*ZT));
	  ZWORK = bench_malloc(N * sizeof(*ZWORK));
	  XT = XWORK = 0;
	  FJCINI(&N, &IO, &IBF, IP, ZT);
     } else {
	  ZT = ZWORK = 0;
	  XT = bench_malloc(3 * N * sizeof(*XT));
	  XWORK = bench_malloc(N * sizeof(*XWORK));
	  FJRINI(&N, &IO, &IBF, IP, XT);
     }
}

void doit(int iter, struct problem *p)
{
     int i;

     if (p->kind == PROBLEM_COMPLEX) {
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
     } else {
	  bench_real *in = p->in;
	  bench_real *out = p->out;
	  if (in == out) {
	       for (i = 0; i < iter; ++i) {
		    FJRRUN(in, XWORK, XWORK, XT, IP);
	       }
	  } else {
	       for (i = 0; i < iter; ++i) {
		    FJRRUN(in, out, XWORK, XT, IP);
	       }
	  }
     }
}

void done(struct problem *p)
{
     if (ZWORK) bench_free(ZWORK);
     if (ZT) bench_free(ZT);
     if (XWORK) bench_free(XWORK);
     if (XT) bench_free(XT);
     UNUSED(p);
}
