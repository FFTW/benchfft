/* this program is in the public domain */

#include "bench-user.h"
#include <stdio.h>
#include <math.h>

#include "cparam.h"

static const char *mknote(void)
{
     static char buf[512];

     sprintf(buf, "Benchmark uses default blocking parameter NBLK=%d, "
	     "padding parameter NP=%d, and cache parameter L2SIZE=%d",
	     NBLK, NP, L2SIZE);
     return buf;
}

BEGIN_BENCH_DOC
BENCH_DOC("name", NAME)
BENCH_DOC("version", "6.0")
BENCH_DOC("author", "Daisuke Takahashi")
BENCH_DOC("year", "2014")
BENCH_DOC("email", "daisuke@is.tsukuba.ac.jp")
BENCH_DOC("language", "Fortran 77")
BENCH_DOC("url", "http://www.ffte.jp/")
BENCH_DOC("url-was-valid-on", "Fri Jul  5 16:41:06 EDT 2019")
DOC_ITEMS
BENCH_DOC("notes", "The backwards transform is scaled.")
BENCH_DOCF("notes", mknote)
BENCH_DOC("bibitem", 
	  "Daisuke Takahashi, A Parallel 3-D FFT Algorithm on Clusters of Vector SMPs, Proc. 5th International Workshop on Applied Parallel Computing (PARA 2000), Lecture Notes in Computer Science, No. 1947, Springer-Verlag, pp. 316-323 (2001).")
BENCH_DOC("bibitem", 
	  "Daisuke Takahashi, A Blocking Algorithm for FFT on Cache-Based Processors, Proc. 9th International Conference on High Performance Computing and Networking Europe (HPCN Europe 2001), Lecture Notes in Computer Science, No. 2110, Springer-Verlag, pp. 551-554 (2001).")
BENCH_DOC("bibitem", 
	  "Daisuke Takahashi, A Blocking Algorithm for Parallel 1-D FFT on Shared-Memory Parallel Computers, Proc. 6th International Conference on Applied Parallel Computing (PARA 2002), Lecture Notes in Computer Science, Springer-Verlag, in press.")
BENCH_DOC("copyright",
	  "Copyright(C) 2000-2002 Daisuke Takahashi\n"
	  "Institute of Information Sciences and Electronics,\n"
	  "University of Tsukuba\n"
	  "1-1-1 Tennodai, Tsukuba-shi, Ibaraki 305-8573, Japan\n"
	  "e-mail: daisuke@is.tsukuba.ac.jp\n"
	  "You may use, copy, modify this code for any purpose and\n"
	  "without fee.  You may distribute this ORIGINAL package.")
END_BENCH_DOC

#define ZFFT1D F77_FUNC(zfft1d,ZFFT1D)
#define ZFFT2D F77_FUNC(zfft2d,ZFFT2D)
#define ZFFT3D F77_FUNC(zfft3d,ZFFT3D)

extern void ZFFT1D(bench_complex *a, 
		   unsigned int *n, int *idir, bench_complex *b);

#ifndef FFTE_VECTOR
#  define FFTE_VECTOR 0
#endif

#if FFTE_VECTOR
extern void ZFFT2D(bench_complex *a, 
		   unsigned int *nx, unsigned int *ny, int *idir,
		   bench_complex *work);
extern void ZFFT3D(bench_complex *a, unsigned int *nx, 
		   unsigned int *ny, unsigned int *nz, int *idir,
		   bench_complex *work);
#  define ZFFT2Dw(a,nx,ny,idir,work) ZFFT2D(a,nx,ny,idir,work)
#  define ZFFT3Dw(a,nx,ny,nz,idir,work) ZFFT3D(a,nx,ny,nz,idir,work)
#else
extern void ZFFT2D(bench_complex *a, 
		   unsigned int *nx, unsigned int *ny, int *idir);
extern void ZFFT3D(bench_complex *a, unsigned int *nx, 
		   unsigned int *ny, unsigned int *nz, int *idir);
#  define ZFFT2Dw(a,nx,ny,idir,work) ZFFT2D(a,nx,ny,idir)
#  define ZFFT3Dw(a,nx,ny,nz,idir,work) ZFFT3D(a,nx,ny,nz,idir)
#endif

int maxdim(const struct problem *p)
{
     unsigned int i, nmax = 0;
     for (i = 0; i < p->rank; ++i)
	  if (p->n[i] > nmax)
	       nmax = p->n[i];
     return nmax;
}

int can_do(struct problem *p)
{
     return (DOUBLE_PRECISION &&
	     p->rank >= 1 && p->rank <= 3 &&
	     p->kind == PROBLEM_COMPLEX &&
	     check_prime_factors(p->size, 5) &&
             problem_in_place(p) &&
	     (p->rank != 2 || maxdim(p) <= NDA2) &&
	     (p->rank != 3 || maxdim(p) <= NDA3)
	  );
}

bench_complex *work = 0;

void setup(struct problem *p)
{
     int idir = 0; // 0 = initialization

     BENCH_ASSERT(can_do(p));
     
     work = (bench_complex *) bench_malloc(sizeof(bench_complex) * p->size
					   * ((p->rank == 1) + FFTE_VECTOR));
     
     /* initialize trig stuff, etc. */
     switch (p->rank) {
	 case 1: 
	      ZFFT1D(p->in, p->n, &idir, work); 
	      break;
	 case 2: 
	      ZFFT2Dw(p->in, p->n + 1, p->n, &idir, work); 
	      break;
	 case 3: 
	      ZFFT3Dw(p->in, p->n + 2, p->n + 1, p->n, &idir, work); 
	      break;
	 default:
	      BENCH_ASSERT(0);
     }
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, 1);
}

void doit(int iter, struct problem *p)
{
     int idir = p->sign;
     int i;

     switch (p->rank) {
	 case 1:
	 {
	      unsigned int *n0 = p->n;
	      for (i = 0; i < iter; ++i)
		   ZFFT1D(p->in, n0, &idir, work); 
	      break;
	 }
	 case 2: 
	 {
	      unsigned int *n0 = p->n;
	      unsigned int *n1 = p->n + 1;
	      for (i = 0; i < iter; ++i)
		   ZFFT2Dw(p->in, n1, n0, &idir, work); 
	      break;
	 }
	 case 3: 
	 {
	      unsigned int *n0 = p->n;
	      unsigned int *n1 = p->n + 1;
	      unsigned int *n2 = p->n + 2;
	      for (i = 0; i < iter; ++i)
		   ZFFT3Dw(p->in, n2, n1, n0, &idir, work); 
	      break;
	 }
	 default:
	      BENCH_ASSERT(0);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(work);
}
