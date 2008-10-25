/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

/* link with -lveclib or -lveclib8 (64-bit) and -lcl ... may also
   need -Wl,-aarchive_shared,-L<PATH> */
/* use MLIB_NUMBER_OF_THREADS to specify parellism (1, or unset, for serial) */

BEGIN_BENCH_DOC
BENCH_DOC("name", "hp-mlib")
BENCH_DOC("package", "HP Mathematical Software Library (mlib, a.k.a. veclib)")
BENCH_DOC("url", "http://h21007.www2.hp.com/portal/site/dspp/menuitem.863c3e4cbcdc3f3515b49c108973a801/?ciid=c008a8ea6ce02110a8ea6ce02110275d6e10RCRD")
BENCH_DOC("url-was-valid-on", "Sat Oct 25 16:20:23 EDT 2008")
BENCH_DOC("notes", "We benchmark the complex (interleaved real/imag) format for complex transforms.")
BENCH_DOC("notes", "We benchmark the non-redundant format for real data (real input array, half of complex output array).")
BENCH_DOC("notes", "We benchmark the unscaled inverse in 1d.")
END_BENCH_DOC

/* defines complex8_t and complex16_t and MLIB_INT: */
#if defined(HAVE_VECLIB8)
#  include <veclib8.h>
#else
#  include <veclib.h>
#endif

/* prefixes for precision */
#ifdef BENCHFFT_SINGLE
#  define C(n) c ## n
#  define S(n) s ## n
typedef complex8_t cmplx;
#else
#  define C(n) z ## n
#  define S(n) d ## n
typedef complex16_t cmplx;
#endif

int can_do(struct problem *p)
{
     return (p->rank >= 1 && p->rank <= 3 &&
	     check_prime_factors(p->size, 5) &&
	     problem_in_place(p) &&
	     /* last dimension must be even for real transforms */
	     (p->kind == PROBLEM_COMPLEX || !(p->n[p->rank - 1] & 1)));
}

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
     copy_r2c_unpacked(p, out);	  
}

void copy_c2r(struct problem *p, bench_complex *in)
{
     copy_c2r_unpacked(p, in);
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     if (p->rank > 1)
	  unnormalize(p, out, 1);
}

bench_real *work = 0;

void setup(struct problem *p)
{
     MLIB_INT iopt = -3, ierr = 0;

     BENCH_ASSERT(can_do(p));
 
     if (p->kind == PROBLEM_COMPLEX) {
	  if (p->rank == 1) {
	       MLIB_INT n1 = p->n[0];
	       work = (bench_real *) bench_malloc(sizeof(bench_real) *
						  5 * (n1/2+1));
	       C(1dfft)((cmplx *) p->in, &n1, work, &iopt, &ierr);
	  }
     } else {
	  if (p->rank == 1) {
	       MLIB_INT n1 = p->n[0];
	       work = (bench_real *) bench_malloc(sizeof(bench_real) * 2 * n1);
	       S(rc1ft)((bench_real *) p->in, &n1, work, &iopt, &ierr);
	  }
     }

     BENCH_ASSERT(ierr == 0);
}

void doit(int iter, struct problem *p)
{
     int i;
     void *in = p->in;
     MLIB_INT iopt = p->sign < 0 ? +1 : -2;
     MLIB_INT ierr;

     if (p->kind == PROBLEM_COMPLEX) {
	  if (p->rank == 1) {
	       MLIB_INT n1 = p->n[0];
               for (i = 0; i < iter; ++i)
		    C(1dfft)((cmplx *) in, &n1, work, &iopt, &ierr);
	  }
	  else if (p->rank == 2) {
	       MLIB_INT n1 = p->n[1];
	       MLIB_INT n2 = p->n[0];
               for (i = 0; i < iter; ++i)
		    C(2dfft)((cmplx *) p->in, &n1, &n2, &n1, &iopt, &ierr);
	  }
	  else if (p->rank == 3) {
	       MLIB_INT n1 = p->n[2];
	       MLIB_INT n2 = p->n[1];
	       MLIB_INT n3 = p->n[0];
               for (i = 0; i < iter; ++i)
		    C(3dfft)((cmplx *) p->in, 
			     &n1,&n2,&n3, &n1,&n2, &iopt, &ierr);
	  }
     } else { /* PROBLEM_REAL */
	  if (p->rank == 1) {
	       MLIB_INT n1 = p->n[0];
               for (i = 0; i < iter; ++i)
		    S(rc1ft)((bench_real *) p->in, &n1, work, &iopt, &ierr);
	  }
	  else if (p->rank == 2) {
	       MLIB_INT n1 = p->n[1];
	       MLIB_INT n2 = p->n[0];
	       MLIB_INT ldx = n1 + 2;
               for (i = 0; i < iter; ++i)
		    S(rc2ft)((bench_real *) p->in,
			     &n1,&n2, &ldx, &iopt, &ierr);
	  }
	  else if (p->rank == 3) {
	       MLIB_INT n1 = p->n[2];
	       MLIB_INT n2 = p->n[1];
	       MLIB_INT n3 = p->n[0];
	       MLIB_INT ldx = n1 + 2;
               for (i = 0; i < iter; ++i)
		    S(rc3ft)((bench_real *) p->in,
			     &n1,&n2,&n3, &ldx,&n2, &iopt, &ierr);
	  }
     }
}

void done(struct problem *p) 
{ 
     UNUSED(p); 
     bench_free(work);
} 
