/* not worth copyrighting */
/* $Id: accopy-to.c,v 1.2 2002-08-15 14:23:58 athena Exp $ */
#include "bench.h"

/* default routine, can be overridden by user */
void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     UNUSED(p);
     UNUSED(out);
}
