/* not worth copyrighting */
/* $Id: accopy-to.c,v 1.1 2001-07-08 20:22:34 athena Exp $ */
#include "bench.h"

/* default routine, can be overridden by user */
void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     UNUSED(p);
     UNUSED(out);
}
