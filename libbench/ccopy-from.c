/* not worth copyrighting */
/* $Id: ccopy-from.c,v 1.3 2001-07-07 22:34:23 athena Exp $ */
#include "bench.h"

/* default routine, can be overridden by user */
void problem_ccopy_from(struct problem *p, bench_complex *in)
{
     if (p->kind == PROBLEM_COMPLEX)
	  cacopy(in, p->in, p->size);
     else {
	  if (p->sign == -1) {
	       /* forward real->hermitian transform */
	       copy_c2r(p, in);
	  } else {
	       /* backward hermitian->real transform */
	       copy_c2h(p, in);
	  }
     }
}
