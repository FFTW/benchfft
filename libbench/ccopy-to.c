/* not worth copyrighting */
/* $Id: ccopy-to.c,v 1.2 2001-07-07 21:56:10 athena Exp $ */
#include "bench.h"

/* default routine, can be overridden by user */
void problem_ccopy_to(struct problem *p, bench_complex *out)
{
     if (p->kind == PROBLEM_COMPLEX)
	  cacopy(p->out, out, p->size);
     else {
	  /* 
	   * user must provide routine because there is no ``standard''
	   * way of doing real transforms
	   */
	  BENCH_ASSERT(((void)"see source code", 0));
     }
}
