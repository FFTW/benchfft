/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "imsl")
BENCH_DOC("package", "IMSL (International Mathematical and Statistical Library)")
BENCH_DOC("author", "Visual Numerics, Inc.")
BENCH_DOC("url", "http://www.vni.com/products/imsl/")
BENCH_DOC("url-was-valid-on", "Mon Jan 14 21:15:19 EST 2002")
BENCH_DOC("notes", "IMSL's FFT routines are based on FFTPACK.")
BENCH_DOC("notes", "Benchmark uses the Fortran 77 interface.")
BENCH_DOC("copyright", "Proprietary, commercial software.")
END_BENCH_DOC

#if defined(BENCHFFT_SINGLE)
#  define F2TRF F77_FUNC(f2trf, F2TRF)
#  define F2TRB F77_FUNC(f2trb, F2TRB)
#  define FFTRI F77_FUNC(fftri, FFTRI)
#  define F2TCF F77_FUNC(f2tcf, F2TCF)
#  define F2TCB F77_FUNC(f2tcb, F2TCB)
#  define FFTCI F77_FUNC(fftci, FFTCI)
#  define F2T2F F77_FUNC(f2t2d, F2T2D)
#  define F2T2B F77_FUNC(f2t2b, F2T2B)
#  define F2T3F F77_FUNC(f2t3f, F2T3F)
#  define F2T3B F77_FUNC(f2t3b, F2T3B)
#else /* double precision */
#  define F2TRF F77_FUNC(df2trf, DF2TRF)
#  define F2TRB F77_FUNC(df2trb, DF2TRB)
#  define FFTRI F77_FUNC(dfftri, DFFTRI)
#  define F2TCF F77_FUNC(df2tcf, DF2TCF)
#  define F2TCB F77_FUNC(df2tcb, DF2TCB)
#  define FFTCI F77_FUNC(dfftci, DFFTCI)
#  define F2T2F F77_FUNC(df2t2d, DF2T2D)
#  define F2T2B F77_FUNC(df2t2b, DF2T2B)
#  define F2T3F F77_FUNC(df2t3f, DF2T3F)
#  define F2T3B F77_FUNC(df2t3b, DF2T3B)
#endif

extern void F2TRF(unsigned int *n, bench_real *x, bench_real *y, bench_real *trig);
extern void F2TRB(unsigned int *n, bench_real *x, bench_real *y, bench_real *trig);
extern void FFTRI(unsigned int *n, bench_real *trig);

extern void F2TCF(unsigned int *n, bench_complex *x, bench_complex *y, bench_real *trig, bench_real *work);
extern void F2TCB(unsigned int *n, bench_complex *x, bench_complex *y, bench_real *trig, bench_real *work);
extern void FFTCI(unsigned int *n, bench_real *trig);

extern void F2T2F(unsigned int *n1, unsigned int *n2, bench_complex *x,  unsigned int *ldx, bench_complex *y, unsigned int *ldy, bench_real *trig1, bench_real *trig2, bench_complex *cwork, bench_real *work);
extern void F2T2B(unsigned int *n1, unsigned int *n2, bench_complex *x,  unsigned int *ldx, bench_complex *y, unsigned int *ldy, bench_real *trig1, bench_real *trig2, bench_complex *cwork, bench_real *work);

extern void F2T3F(unsigned int *n1, unsigned int *n2, unsigned int *n3, bench_complex *x,  unsigned int *ldx, unsigned int *mdx, bench_complex *y, unsigned int *ldy, unsigned int *mdy, bench_real *trig1, bench_real *trig2, bench_real *trig3, bench_real *work);
extern void F2T3B(unsigned int *n1, unsigned int *n2, unsigned int *n3, bench_complex *x,  unsigned int *ldx, unsigned int *mdx, bench_complex *y, unsigned int *ldy, unsigned int *mdy, bench_real *trig1, bench_real *trig2, bench_real *trig3, bench_real *work);

int can_do(struct problem *p)
{
     return (p->rank >= 1 &&
	     ((p->rank <= 3 && p->kind == PROBLEM_COMPLEX) || (p->rank == 1)));
}

bench_real *trig = 0, *trig2 = 0, *trig3 = 0, *work = 0;

#define MAX2(a,b) ((a) > (b) ? (a) : (b))
#define MAX3(a,b,c) MAX2(a,MAX2(b,c))

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));

     if (p->kind == PROBLEM_COMPLEX) {
	  if (p->rank == 1) {
	       trig = (bench_real*) bench_malloc(sizeof(bench_real) * 
						 (4*p->size + 15));
	       FFTCI(&p->size, trig);

	       work = (bench_real*) bench_malloc(sizeof(bench_real) *
						 2 * p->size);
	  }
	  else if (p->rank == 2) {
	       trig = (bench_real*) bench_malloc(sizeof(bench_real) *
                                                 (4*p->n[1] + 15));
               FFTCI(&p->n[1], trig);

	       trig2 = (bench_real*) bench_malloc(sizeof(bench_real) *
						  (4*p->n[0] + 15));
               FFTCI(&p->n[0], trig2);

	       work = (bench_real*) bench_malloc(sizeof(bench_real) *
						 2 * MAX2(p->n[0],
							  p->n[1]));
	  }
	  else /* (p->rank == 3) */ {
	       trig = (bench_real*) bench_malloc(sizeof(bench_real) *
                                                 (4*p->n[2] + 15));
               FFTCI(&p->n[2], trig);

	       trig2 = (bench_real*) bench_malloc(sizeof(bench_real) *
						  (4*p->n[1] + 15));
               FFTCI(&p->n[1], trig2);

	       trig3 = (bench_real*) bench_malloc(sizeof(bench_real) *
						  (4*p->n[0] + 15));
               FFTCI(&p->n[0], trig3);

	       work = (bench_real*) bench_malloc(sizeof(bench_real) *
						 2 * MAX3(p->n[0],
							  p->n[1],
							  p->n[2]));
	  }
     }
     else /* PROBLEM_REAL */ {
	  trig = (bench_real*) bench_malloc(sizeof(bench_real) * 
					    (2*p->size + 15));
	  FFTRI(&p->size, trig);
     }
}

void doit(int iter, struct problem *p)
{
     int i;

     if (p->kind == PROBLEM_COMPLEX) {
	  bench_complex *x = (bench_complex *) p->in;
	  bench_complex *y = (bench_complex *) p->out;
	  switch (p->rank) {
	      case 1:
	      {
		   unsigned int n = p->size;
		   if (p->sign == -1)
			for (i = 0; i < iter; ++i)
			     F2TCF(&n, x, y, trig, work);
		   else
			for (i = 0; i < iter; ++i)
			     F2TCB(&n, x, y, trig, work);
		   break;
	      }
	      case 2:
	      {
		   unsigned int n1 = p->n[1], n2 = p->n[0];
		   bench_complex cwork;
		   if (p->sign == -1)
			for (i = 0; i < iter; ++i)
			     F2T2F(&n1, &n2, x, &n1, y, &n1,
				   trig, trig2, &cwork, work);
		   else
			for (i = 0; i < iter; ++i)
			     F2T2B(&n1, &n2, x, &n1, y, &n1,
				   trig, trig2, &cwork, work);
		   break;
	      }
	      case 3:
	      {
		   unsigned int n1 = p->n[2], n2 = p->n[1], n3 = p->n[0];
		   if (p->sign == -1)
			for (i = 0; i < iter; ++i)
			     F2T3F(&n1, &n2, &n3, x, &n1, &n2, y, &n1, &n2,
				   trig, trig2, trig3, work);
		   else
			for (i = 0; i < iter; ++i)
			     F2T3B(&n1, &n2, &n3, x, &n1, &n2, y, &n1, &n2,
				   trig, trig2, trig3, work);
		   break;
	      }
	      default:
		   BENCH_ASSERT(0);
	  }
     }
     else /* PROBLEM_REAL */ {
	  bench_real *x = (bench_real *) p->in;
	  bench_real *y = (bench_real *) p->out;
	  unsigned int n = p->size;
	  if (p->sign == -1)
	       for (i = 0; i < iter; ++i)
		    F2TRF(&n, x, y, trig);
	  else
	       for (i = 0; i < iter; ++i)
		    F2TRB(&n, x, y, trig);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(work);
     bench_free(trig3);
     bench_free(trig2);
     bench_free(trig);
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_1d_fftpack(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_fftpack(p, in, -1.0);
}
