/* not worth copyrighting */
/* $Id: ccopy-to.c,v 1.1 2001-07-07 14:16:56 athena Exp $ */
#include "bench.h"

/* default routine, can be overridden by user */
void problem_ccopy_to(struct problem *p, bench_complex *out)
{
     cacopy(p->out, out, p->size);
}
