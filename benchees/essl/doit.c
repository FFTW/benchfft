/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "essl")
BENCH_DOC("author", "IBM")
BENCH_DOC("copyright", "Proprietary software.")
BENCH_DOC("notes", "ESSL is the IBM Engineering and Scientific Subroutine Library.")
BENCH_DOC("url", "http://www.rs6000.ibm.com/software/Apps/essl.html")
BENCH_DOC("url", "http://www.rs6000.ibm.com/doc_link/en_US/a_doc_lib/sp32/essl/essl065.html")
BENCH_DOC("url-was-valid-on", "Tue Aug 14 15:07:49 EDT 2001")
END_BENCH_DOC

#ifdef BENCHFFT_SINGLE
#  define TYPE_PREFIX s
#else
#  define TYPE_PREFIX d
#endif

#define PASTEx(a,b) a##b
#define PASTE(a,b) PASTEx(a,b) /* need this hack to paste symbols */

#define CFT PASTE(TYPE_PREFIX, cft)
extern void CFT(int *init,
		bench_complex *x, int *xstride, int *xdist,
		bench_complex *y, int *ystride, int *ydist,
		int *n, int *howmany, int *isign, bench_real *scale,
		double *work1, int *nwork1, 
		double *work2, int *nwork2);

#define CFT2 PASTE(TYPE_PREFIX, cft2)
extern void CFT2(int *init,
		 bench_complex *x, int *xstride1, int *xstride2,
		 bench_complex *y, int *ystride1, int *ystride2,
		 int *n1, int *n2, int *isign, bench_real *scale,
		 double *work1, int *nwork1, 
		 double *work2, int *nwork2);

#define CFT3 PASTE(TYPE_PREFIX, cft3)
extern void CFT3(bench_complex *x, int *xstride2, int *xstride3,
		 bench_complex *y, int *ystride2, int *ystride3,
		 int *n1, int *n2, int *n3, int *isign, bench_real *scale,
		 double *work, int *nwork);

#define RCFT PASTE(TYPE_PREFIX, rcft)
extern void RCFT(int *init,
		 bench_real *x, int *xdist,
		 bench_complex *y, int *ydist,
		 int *n, int *howmany, int *isign, bench_real *scale,
		 double *work1, int *nwork1, 
		 double *work2, int *nwork2);

#define CRFT PASTE(TYPE_PREFIX, crft)
extern void CRFT(int *init,
		 bench_complex *y, int *ydist,
		 bench_real *x, int *xdist,
		 int *n, int *howmany, int *isign, bench_real *scale,
		 double *work1, int *nwork1, 
		 double *work2, int *nwork2);

#define RCFT2 PASTE(TYPE_PREFIX, rcft2)
extern void RCFT2(int *init,
		  bench_real *x, int *xstride2,
		  bench_complex *y, int *ystride2,
		  int *n1, int *n2, int *isign, bench_real *scale,
		  double *work1, int *nwork1, 
		  double *work2, int *nwork2);

#define RCFT3 PASTE(TYPE_PREFIX, rcft3)
extern void RCFT3(bench_real *x, int *xstride2, int *xstride3,
		  bench_complex *y, int *ystride2, int *ystride3,
		  int *n1, int *n2, int *n3, int *isign, bench_real *scale,
		  double *work, int *nwork);

#define CRFT2 PASTE(TYPE_PREFIX, crft2)
extern void CRFT2(int *init,
		  bench_complex *y, int *ystride2,
		  bench_real *x, int *xstride2,
		  int *n1, int *n2, int *isign, bench_real *scale,
		  double *work1, int *nwork1, 
		  double *work2, int *nwork2);

#define CRFT3 PASTE(TYPE_PREFIX, crft3)
extern void CRFT3(bench_complex *y, int *ystride2, int *ystride3,
		  bench_real *x, int *xstride2, int *xstride3,
		  int *n1, int *n2, int *n3, int *isign, bench_real *scale,
		  double *work, int *nwork);

int n_ok(int rank, unsigned int *n, int is_real)
{
     int i;
     for (i = 0; i < rank; ++i) {
	  int log2;
	  unsigned int N = n[i];

	  if (N > 37748736)
	       return 0;

	  for (log2 = 0; N > 1 && N % 2 == 0; N /= 2, log2 += 1);
	  if (log2 < 1 || log2 > 25)
	       return 0;

	  if (!check_prime_factors(N, 11))
	       return 0;

	  if (N % (3*3*3) == 0 ||
	      N % (5*5) == 0 || N % (7*7) == 0 || N % (11*11) == 0)
	       return 0;
     }
     return 1;
}

int can_do(struct problem *p)
{
     return (p->rank >= 1 && p->rank <= 3 && 
	     n_ok(p->rank, p->n, p->kind == PROBLEM_REAL));
}

static int imax2(int a, int b) { return (a > b ? a : b); }
static int imin2(int a, int b) { return (a < b ? a : b); }

double *work1, *work2;
int nwork1, nwork2;

#define MAX2(a,b) ((a) > (b) ? (a) : (b))
#define MIN2(a,b) ((a) < (b) ? (a) : (b))

void setup(struct problem *p)
{
     int init = 1, stride1 = 1, stride2 = p->size, stride1o, stride2o;
     int isign = (p->sign < 0) ? 1 : -1;
     int n1, n2 = 1, n3, nmax;
     bench_real scale = 1.0;

     /* The ESSL authors, in their wisdom, provided a cornucopia of
        lovely formulas for workspace sizes, in lieu of a routine to
        return the correct size. */
     if (p->kind == PROBLEM_COMPLEX) {
	  switch (p->rank) {
	      case 1:
		   n1 = p->n[0];
		   if (SINGLE_PRECISION) {
			if (p->n[0] <= 8192)
			     nwork1 = nwork2 = 20000;
			else
			     nwork1 = nwork2 = 20000 + 1.14 * n1;
		   }
		   else /* (DOUBLE_PRECISION) */ {
			if (p->n[0] <= 2048)
			     nwork1 = nwork2 = 20000;
			else
			     nwork1 = nwork2 = 20000 + 2.28 * n1;
		   }
		   break;

	      case 2:
		   n2 = p->n[0];
		   n1 = p->n[1];
		   stride1 = 1; stride2 = n1;
		   nmax = MAX2(n1,n2);
		   if (SINGLE_PRECISION) {
			if (nmax <= 8192)
			     nwork1 = 40000;
			else
			     nwork1 = 40000 + 1.14 * (n1 + n2);
			if (nmax < 252)
			     nwork2 = 20000;
			else
			     nwork2 = 20000 +
				  (nmax + 256) * (MIN2(MIN2(n1,n2),64)+1.14);
		   }
		   else /* (DOUBLE_PRECISION) */ {
			if (nmax <= 2048)
			     nwork1 = 40000;
			else
			     nwork1 = 40000 + 2.28 * (n1 + n2);
			if (nmax < 252)
			     nwork2 = 20000;
			else
			     nwork2 = 20000 +
				  (2*nmax + 256) * (MIN2(MIN2(n1,n2),64)+1.14);
		   }
		   break;

	      case 3:
		   n3 = p->n[0];
		   n2 = p->n[1];
		   n1 = p->n[2];
		   stride1 = n1; stride2 = n1*n2;
		   nmax = MAX2(n2,n3);
		   nwork2 = 1;
		   if (SINGLE_PRECISION) {
			if (n1 <= 8192)
			     nwork1 = 60000;
			else
			     nwork1 = 60000 + 2.28*n1;
			if (nmax >= 252)
			     nwork1 += (nmax + 256)
				  * (MIN2(64,n1 * (n3 < 252 ? 1 : n2)) + 2.28);
		   }
		   else /* (DOUBLE_PRECISION) */ {
			if (n1 <= 2048)
			     nwork1 = 60000;
			else
			     nwork1 = 60000 + 4.56*n1;
			if (nmax >= 252)
			     nwork1 += (2*nmax + 256)
				  * (MIN2(64,n1 * (n3 < 252 ? 1 : n2)) + 4.56);
		   }
		   break;
	  }
     }
     else if (p->kind == PROBLEM_REAL) {
	  switch (p->rank) {
	      case 1:
		   n1 = p->n[0];
		   stride1o = 1; stride2o = n1/2 + 1;
		   if (SINGLE_PRECISION) {
			if (p->n[0] <= 16384) {
			     nwork1 = 25000;
			     nwork2 = 20000;
			}
			else {
			     nwork1 = 20000 + 0.82 * n1;
			     nwork2 = 20000 + 0.57 * n1;
			}
		   }
		   else /* (DOUBLE_PRECISION) */ {
			if (p->n[0] <= 4096) {
			     nwork1 = 22000;
			     nwork2 = 20000;
			}
			else {
			     nwork1 = 20000 + 1.64 * n1;
			     nwork2 = 20000 + 1.14 * n1;
			}
		   }
		   break;

	      case 2:
		   n2 = p->n[0];
		   n1 = p->n[1];
		   stride1 = 1; stride2 = n1;
		   stride1o = 1; stride2o = n1/2 + 1;
		   nmax = MAX2(n1/2,n2);
		   if (SINGLE_PRECISION) {
			if (nmax <= 8192)
			     nwork1 = 45000;
			else
			     nwork1 = 40000 + 0.82*n1 + 1.14*n2;
			nwork2 = 20000;
			if (n1 > 16384)
			     nwork2 += 0.57*n1;
			if (n2 >= 252)
			     nwork2 += (n2 + 256) * (1.14 + MIN2(64,1+n1/2));
		   }
		   else /* (DOUBLE_PRECISION) */ {
			if (nmax <= 2048)
			     nwork1 = 42000;
			else
			     nwork1 = 40000 + 1.64*n1 + 2.28*n2;
			nwork2 = 20000;
			if (n1 > 4096)
			     nwork2 += 1.14*n1;
			if (n2 >= 252)
			     nwork2 += (2*n2 + 256) * (2.28 + MIN2(64,1+n1/2));
		   }
		   break;

	      case 3:
		   n3 = p->n[0];
		   n2 = p->n[1];
		   n1 = p->n[2];
		   stride1 = n1; stride2 = n1*n2;
		   stride1o = n1/2 + 1; stride2o = stride1o * n2;
		   nmax = MAX2(n2,n3);
		   nwork2 = 1;
		   if (SINGLE_PRECISION) {
			if (n1 <= 16384)
			     nwork1 = 65000;
			else
			     nwork1 = 60000 + 1.39*n1;
			if (nmax >= 252)
			     nwork1 += (nmax + 256)
				  * (MIN2(64, stride1o 
					  * (n3 < 252 ? 1 : n2)) + 2.28);
		   }
		   else /* (DOUBLE_PRECISION) */ {
			if (n1 <= 4096)
			     nwork1 = 62000;
			else
			     nwork1 = 60000 + 2.78*n1;
			if (nmax >= 252)
			     nwork1 += (2*nmax + 256)
				  * (MIN2(64, stride1o 
					  * (n3 < 252 ? 1 : n2)) + 4.56);
		   }
		   break;
	  }
	  if (problem_in_place(p)) {
	       stride1 = stride1o == 1 ? 1 : stride1o * 2;
	       stride2 = stride2o * 2;
	  }
     }

     work1 = (double *) bench_malloc(sizeof(double) * nwork1);
     work2 = (double *) bench_malloc(sizeof(double) * nwork2);

     if (p->kind == PROBLEM_COMPLEX) {
	  switch (p->rank) {
	      case 1:
		   CFT(&init,
		       (bench_complex *) p->in, &stride1, &stride2, 
		       (bench_complex *) p->out, &stride1, &stride2,
		       &n1, &n2, &isign, &scale,
		       work1, &nwork1, work2, &nwork2);
		   break;
	      case 2:
		   CFT2(&init,
			(bench_complex *) p->in, &stride1, &stride2, 
			(bench_complex *) p->out, &stride1, &stride2,
			&n1, &n2, &isign, &scale,
			work1, &nwork1, work2, &nwork2);
		   break;
	      case 3:
		   break;
	  }
     }
     else if (p->kind == PROBLEM_REAL) {
	  switch (p->rank) {
	      case 1:
		   if (p->sign == -1)
			RCFT(&init,
			     (bench_real *) p->in, &stride2, 
			     (bench_complex *) p->out, &stride2o,
			     &n1, &n2, &isign, &scale,
			     work1, &nwork1, work2, &nwork2);
		   else
			CRFT(&init,
			     (bench_complex *) p->in, &stride2o, 
			     (bench_real *) p->out, &stride2,
			     &n1, &n2, &isign, &scale,
			     work1, &nwork1, work2, &nwork2);
		   break;
	      case 2:
		   if (p->sign == -1)
			RCFT2(&init,
			      (bench_real *) p->in, &stride2, 
			      (bench_complex *) p->out, &stride2o,
			      &n1, &n2, &isign, &scale,
			      work1, &nwork1, work2, &nwork2);
		   else
			CRFT2(&init,
			      (bench_complex *) p->in, &stride2o,
			      (bench_real *) p->out, &stride2, 
			      &n1, &n2, &isign, &scale,
			      work1, &nwork1, work2, &nwork2);
		   break;
	      case 3:
		   break;
	  }
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     int init = 0;
     int isign = (p->sign < 0) ? 1 : -1;
     int n1, n2, n3, stride1, stride2, stride1o, stride2o;
     bench_real scale = 1.0;

      if (p->kind == PROBLEM_COMPLEX) {
	  switch (p->rank) {
	      case 1:
		   n1 = p->n[0]; n2 = 1;
		   stride1 = 1; stride2 = n1;
		   for (i = 0; i < iter; ++i) {
			CFT(&init,
			    (bench_complex *) p->in, &stride1, &stride2, 
			    (bench_complex *) p->out, &stride1, &stride2,
			    &n1, &n2, &isign, &scale,
			    work1, &nwork1, work2, &nwork2);
		   }
		   break;
	      case 2:
		   n1 = p->n[1]; n2 = p->n[0];
		   stride1 = 1; stride2 = n1;
		   for (i = 0; i < iter; ++i) {
			CFT2(&init,
			     (bench_complex *) p->in, &stride1, &stride2, 
			     (bench_complex *) p->out, &stride1, &stride2,
			     &n1, &n2, &isign, &scale,
			     work1, &nwork1, work2, &nwork2);
		   }
		   break;
	      case 3:
		   n1 = p->n[2]; n2 = p->n[1]; n3 = p->n[0];
		   stride1 = n1; stride2 = n1*n2;
		   for (i = 0; i < iter; ++i) {
			CFT3((bench_complex *) p->in, &stride1, &stride2, 
			     (bench_complex *) p->out, &stride1, &stride2,
			     &n1, &n2, &n3, &isign, &scale,
			     work1, &nwork1);
		   }
		   break;
	  }
     }
      else if (p->kind == PROBLEM_REAL) {
	  switch (p->rank) {
	      case 1:
		   n1 = p->n[0]; n2 = 1;
		   stride2 = n1; stride2o = n1/2 + 1;
		   if (problem_in_place(p))
			stride2 = stride2o * 2;
		   if (p->sign == -1)
			for (i = 0; i < iter; ++i) {
			     RCFT(&init,
				  (bench_real *) p->in, &stride2, 
				  (bench_complex *) p->out, &stride2,
				  &n1, &n2, &isign, &scale,
				  work1, &nwork1, work2, &nwork2);
			}
		   else
			for (i = 0; i < iter; ++i) {
			     CRFT(&init,
				  (bench_complex *) p->in, &stride2,
				  (bench_real *) p->out, &stride2, 
				  &n1, &n2, &isign, &scale,
				  work1, &nwork1, work2, &nwork2);
			}
		   break;
	      case 2:
		   n1 = p->n[1]; n2 = p->n[0];
		   stride1 = 1; stride2 = n1;
		   stride2o = n1/2 + 1;
		   if (problem_in_place(p))
                        stride2 = stride2o * 2;
		   if (p->sign == -1)
			for (i = 0; i < iter; ++i) {
			     RCFT2(&init,
				   (bench_real *) p->in, &stride2, 
				   (bench_complex *) p->out, &stride2o,
				   &n1, &n2, &isign, &scale,
				   work1, &nwork1, work2, &nwork2);
			}
		   else
			for (i = 0; i < iter; ++i) {
			     CRFT2(&init,
				   (bench_complex *) p->in, &stride2o,
				   (bench_real *) p->out, &stride2, 
				   &n1, &n2, &isign, &scale,
			      work1, &nwork1, work2, &nwork2);
			}
		   break;
	      case 3:
		   n1 = p->n[2]; n2 = p->n[1]; n3 = p->n[0];
		   stride1 = n1; stride2 = n1*n2;
		   stride1o = n1/2 + 1; stride2o = stride1o * n2;
		   if (problem_in_place(p)) {
			stride1 = stride1o * 2;
                        stride2 = stride2o * 2;
		   }
		   if (p->sign == -1)
			for (i = 0; i < iter; ++i) {
			     RCFT3((bench_real *) p->in, &stride1, &stride2, 
				   (bench_complex *) p->out,
				   &stride1o, &stride2o,
				   &n1, &n2, &n3, &isign, &scale,
				   work1, &nwork1);
			}
		   else
			for (i = 0; i < iter; ++i) {
			     CRFT3((bench_complex *) p->in,
				   &stride1o, &stride2o,
				   (bench_real *) p->out, &stride1, &stride2, 
				   &n1, &n2, &n3, &isign, &scale,
				   work1, &nwork1);
			}
		   break;
	  }
     }
}

void done(struct problem *p)
{
     bench_free(work2);
     bench_free(work1);
}

/* The ESSL real transforms use the same format as rfftwnd. */

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_unpacked(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_unpacked(p, in, -1.0);
}

void copy_r2c(struct problem *p, bench_complex *out)
{
     if (problem_in_place(p))
	  copy_r2c_unpacked(p, out);	  
     else
	  copy_r2c_packed(p, out);
}

void copy_c2r(struct problem *p, bench_complex *in)
{
     if (problem_in_place(p))
	  copy_c2r_unpacked(p, in);
     else
	  copy_c2r_packed(p, in);
}
