/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "cwplib")
BENCH_DOC("author", "Dave Hale")
BENCH_DOC("year", "1989")
BENCH_DOC("language", "C")
BENCH_DOC("url", "http://risc2.numis.nwu.edu/ftp/pub/transforms/cwplib.tar.gz")
BENCH_DOC("url-was-valid-on", "Sat Aug 31 16:39:12 EDT 2002")
BENCH_DOC("copyright", "Copyright (c) Colorado School of Mines, 1995.")
END_BENCH_DOC

#include "cwp.h"

int can_do(struct problem *p)
{
     unsigned int i;

     if (!SINGLE_PRECISION)
	  return 0;

     if (p->kind == PROBLEM_COMPLEX) {
	  /* 
	   * CWPLIB does not really provide a 2D transform.
	   * Nevertheless, source code comments provide instructions
	   * on how to use CWPLIB to compute a 2D transform, which we
	   * take as an implementation.
	   */
	  if (p->rank != 1 && p->rank != 2)
	       return 0;

	  /* must be in place */
	  if (!problem_in_place(p))
	       return 0;

	  for (i = 0; i < p->rank; ++i) {
	       if ((unsigned int)npfa(p->n[i]) != p->n[i])
		    return 0;
	  } 
     } else {
	  /* no 2D instructions are provided */
	  if (p->rank != 1)
	       return 0;

	  /* can be either in place or out of place */

	  if ((unsigned int)npfar(p->n[0]) != p->n[0])
	       return 0;
     }

     return 1;
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_unpacked(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_unpacked(p, in, -1.0);
}

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
     /* nothing to do */
}

void doit(int iter, struct problem *p)
{
     int i;
     int sign = p->sign;
     void *in = p->in;

     if (p->kind == PROBLEM_COMPLEX) {
	  switch (p->rank) {
	      case 1:
	      {
		   int n = p->n[0];
		   for (i = 0; i < iter; ++i) 
			pfacc(sign, n, in);
		   break;
	      }
	      case 2:
	      {
		   int n0 = p->n[0];
		   int n1 = p->n[1];
		   for (i = 0; i < iter; ++i) { 
			pfa2cc(sign, 1, n1, n0, in);
			pfa2cc(sign, 2, n1, n0, in);
		   }
		   break;
	      }
     
	      default:
		   BENCH_ASSERT(0);
	  }
     } else {
	  /* p->rank == 1*/
	  int n = p->n[0];
	  void *out = p->out;

	  if (sign == -1) {
	       /* real->complex */
	       for (i = 0; i < iter; ++i) 
		    pfarc(sign, n, in, out);
	  } else {
	       for (i = 0; i < iter; ++i) 
		    pfacr(sign, n, in, out);
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
