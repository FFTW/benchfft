/* not worth copyrighting */

/* $Id: deallocate.c,v 1.2 2002-08-15 14:23:58 athena Exp $ */

#include "config.h"
#include "bench.h"

void problem_free(struct problem *p)
{
     if (p->out && p->out != p->in)
	  bench_free(p->out);
     if (p->in)
	  bench_free(p->in);
}
