/* not worth copyrighting */
/* $Id: ccopy-from.c,v 1.1 2001-07-07 14:16:56 athena Exp $ */
#include "bench.h"

/* default routine, can be overridden by user */
void problem_ccopy_from(struct problem *p, bench_complex *in)
{
     cacopy(in, p->in, p->size);
}
