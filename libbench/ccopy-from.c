/* not worth copyrighting */
/* $Id: ccopy-from.c,v 1.2 2001-07-07 21:56:10 athena Exp $ */
#include "bench.h"

/* default routine, can be overridden by user */
void problem_ccopy_from(struct problem *p, bench_complex *in)
{
     if (p->kind == PROBLEM_COMPLEX)
	  cacopy(in, p->in, p->size);
     else {
	  /* 
	   * user must provide routine because there is no ``standard''
	   * way of doing real transforms
	   */
	  BENCH_ASSERT(((void)"see source code", 0));
     }
}
