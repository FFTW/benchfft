/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>
#include <stdio.h>

static const char *mkvers(void)
{
     int major, minor;
     static char buf[160];
     acmlversion(&major, &minor);
     sprintf(buf, "%d.%d", major, minor);
     return buf;
}


BEGIN_BENCH_DOC
BENCH_DOC("name", "acml")
BENCH_DOC("package", "AMD Core Math Library (ACML)")
BENCH_DOCF("version", mkvers)
BENCH_DOC("notes", "transform is scaled by sqrt(n)")
BENCH_DOC("notes", "backward real transform is scaled by 2")
END_BENCH_DOC


int can_do(struct problem *p)
{
     switch (p->kind) {
	 case PROBLEM_COMPLEX:
	      return (p->rank >= 1 && p->rank <= 3 && problem_in_place(p));
	 case PROBLEM_REAL:
	      return (p->rank == 1 && problem_in_place(p));
     }
     return 0;
}

static void *WORK;

#ifdef BENCHFFT_SINGLE
#define CFFT1D cfft1d_
#define CFFT2D cfft2d_
#define CFFT3D cfft3d_
#define SCFFT scfft_
#define CSFFT csfft_
#else
#define CFFT1D zfft1d_
#define CFFT2D zfft2d_
#define CFFT3D zfft3d_
#define SCFFT dzfft_
#define CSFFT zdfft_
#endif


/* Undo normalization: */
void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     bench_complex x;
     c_re(x) = sqrt(p->size);
     c_im(x) = 0;
     
     cascale(out, p->size, x);
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_1d_halfcomplex(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_halfcomplex(p, in, +1.0);
     ascale(p->in, p->size, 2.0);
}

void setup(struct problem *p)
{
     int mode = 0;
     int info;
 
     BENCH_ASSERT(can_do(p));
 
     if (p->kind == PROBLEM_COMPLEX) {
	  switch (p->rank) {
	      case 1:
		   WORK = bench_malloc((2*p->n[0] + 100) *
				       sizeof(bench_complex));
		   CFFT1D(&mode, p->n, p->in, WORK, &info);
		   break;
	      case 2:
		   WORK = bench_malloc((p->n[0] * p->n[1] 
					+ 2 * p->n[0] 
					+ 2 * p->n[1])
				       * sizeof(bench_complex));
		   break;
	      case 3:
		   WORK = bench_malloc((p->n[0] * p->n[1] * p->n[2]
					+ 2 * p->n[0] 
					+ 2 * p->n[1]
					+ 2 * p->n[2])
				       * sizeof(bench_complex));
		   break;
	  }
     } else {
	  WORK = bench_malloc((2*p->n[0] + 100) * sizeof(bench_complex));
	  if (p->sign == -1) {
	       SCFFT(&mode, p->n, p->in, WORK, &info);
	  } else {
	       CSFFT(&mode, p->n, p->in, WORK, &info);
	  }
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     void *in = p->in;
     int sign = p->sign;
     int info;

     if (p->kind == PROBLEM_COMPLEX) {
	  switch (p->rank) {
	      case 1:
		   for (i = 0; i < iter; ++i) 
			CFFT1D(&sign, p->n, p->in, WORK, &info);
		   break;
	      case 2:
		   for (i = 0; i < iter; ++i) 
			CFFT2D(&sign, p->n + 1, p->n, p->in, WORK, &info);
		   break;
	      case 3:
		   for (i = 0; i < iter; ++i) 
			CFFT3D(&sign, p->n + 2, p->n + 1, p->n, 
			       p->in, WORK, &info);
		   break;
	  }
     } else {
	  if (sign == -1) {
	       sign = 1;
	       for (i = 0; i < iter; ++i) 
		    SCFFT(&sign, p->n, p->in, WORK, &info);
	  } else {
	       for (i = 0; i < iter; ++i) 
		    CSFFT(&sign, p->n, p->in, WORK, &info);
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(WORK);
}
