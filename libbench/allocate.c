/* not worth copyrighting */

/* $Id: allocate.c,v 1.1 2001-07-07 14:16:56 athena Exp $ */

#include "config.h"
#include "bench.h"

/*
 * Allocate I/O arrays for a problem.
 *
 * This is the default routine that can be overridden by the user in
 * complicated cases.
 */
void problem_alloc(struct problem *p, int in_place)
{
     size_t s = p->kind == PROBLEM_COMPLEX ? 
	  sizeof(bench_complex) : sizeof(bench_real);
     
     p->in = bench_malloc(p->size * s);

     if (in_place)
	  p->out = p->in;
     else
	  p->out = bench_malloc(p->size * s);
}
