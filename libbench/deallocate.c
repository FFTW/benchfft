/* not worth copyrighting */

/* $Id: deallocate.c,v 1.3 2002-08-16 12:09:55 athena Exp $ */

#include "config.h"
#include "bench.h"

void problem_free(struct problem *p)
{
     if (p->out && p->out != p->in)
	  bench_free(p->out);
     if (p->in)
	  bench_free(p->in);
}
