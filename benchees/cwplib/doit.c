/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "cwplib")
BENCH_DOC("author", "Dave Hale")
BENCH_DOC("year", "1989")
BENCH_DOC("language", "C")
BENCH_DOC("url", "http://risc1.numis.nwu.edu/ftp/pub/transforms/cwplib.tar.gz")
BENCH_DOC("copyright", "Copyright (c) Colorado School of Mines, 1995.")
END_BENCH_DOC

#include "cwp.h"

int can_do(struct problem *p)
{
     unsigned int i;

     if (!SINGLE_PRECISION)
	  return 0;
     if (p->rank != 1 && p->rank != 2)
	  return 0;

     if (p->kind == PROBLEM_COMPLEX) {
	  /* must be in place */
	  if (p->in != p->out)
	       return 0;

	  for (i = 0; i < p->rank; ++i) {
	       if ((unsigned int)npfa(p->n[i]) != p->n[i])
		    return 0;
	  } 
     } else {
	  for (i = 0; i < p->rank; ++i) {
	       if ((unsigned int)npfar(p->n[i]) != p->n[i])
		    return 0;
	  } 
     }

     return 1;
}

/* override default problem_alloc */
void problem_alloc(struct problem *p, int in_place)
{
     if (p->kind == PROBLEM_COMPLEX) {
	  p->in = bench_malloc(p->size * sizeof(bench_complex));
	  
	  if (in_place)
	       p->out = p->in;
	  else
	       p->out = bench_malloc(p->size * sizeof(bench_complex));
     } else {
	  /* TODO */
     }
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
}

void done(struct problem *p)
{
     UNUSED(p);
}
