/* not worth copyrighting */

/* $Id: allocate.c,v 1.2 2001-07-13 00:02:10 athena Exp $ */

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
     if (p->kind == PROBLEM_COMPLEX) {
	  p->in = bench_malloc(p->size * sizeof(bench_complex));
	  
	  if (in_place)
	       p->out = p->in;
	  else
	       p->out = bench_malloc(p->size * sizeof(bench_complex));
     } else {
	  size_t s = 1;
	  unsigned int i;

	  for (i = 0; i < p->rank; ++i)
	       /* slightly overallocate to account for unpacked formats */
	       s *= p->n[i] + 2;

	  p->in = bench_malloc(s * sizeof(bench_real));
	  
	  if (in_place)
	       p->out = p->in;
	  else
	       p->out = bench_malloc(s * sizeof(bench_real));
     }
}
