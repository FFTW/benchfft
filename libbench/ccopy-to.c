/* not worth copyrighting */
/* $Id: ccopy-to.c,v 1.3 2001-07-07 22:34:23 athena Exp $ */
#include "bench.h"

/* default routine, can be overridden by user */
void problem_ccopy_to(struct problem *p, bench_complex *out)
{
     if (p->kind == PROBLEM_COMPLEX)
	  cacopy(p->out, out, p->size);
     else {
	  if (p->sign == -1) {
	       /* forward real->hermitian transform */
	       copy_h2c(p, out);
	  } else {
	       /* backward hermitian->real transform */
	       copy_r2c(p, out);
	  }
     }
}
