/* not worth copyrighting */

/* $Id: deallocate.c,v 1.1 2001-07-07 14:16:56 athena Exp $ */

#include "config.h"
#include "bench.h"

void problem_free(struct problem *p)
{
     if (p->out && p->out != p->in)
	  bench_free(p->out);
     if (p->in)
	  bench_free(p->in);
}
