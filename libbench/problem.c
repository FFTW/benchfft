/*
 * Copyright (c) 2001 Matteo Frigo
 * Copyright (c) 2001 Steven G. Johnson
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* $Id: problem.c,v 1.1 2001-07-04 22:50:39 athena Exp $ */

#include "config.h"
#include "bench.h"
#include <stdio.h>
#include <stdlib.h>


/* parse a problem description, return a problem */
/* TODO: only knows complex 1D for now */
struct problem *problem_parse(const char *desc)
{
     int n;
     int in_place = 0;
     int sign = -1;
     struct problem *p;
     const bench_complex czero = {0, 0};

     for (;;) {
	  switch (*desc) {
	      case 'i':
	      case 'I':
		   in_place = 1;
		   ++desc;
		   continue;
	      case 'f':
	      case 'F':
		   sign = -1;
		   ++desc;
		   continue;
	      case 'b':
	      case 'B':
		   sign = 1;
		   ++desc;
		   continue;
	      default:
		   ;
	  }
	  break;
     }

     n = atoi(desc);

     p = bench_malloc(sizeof(struct problem));

     p->kind = PROBLEM_COMPLEX;
     p->userinfo = 0;
     p->rank = 1;
     p->n = bench_malloc(p->rank * sizeof(int));
     p->n[0] = n;
     p->size = n;

     p->p.complex.sign = sign;
     p->p.complex.in = bench_malloc(n * sizeof(bench_complex));

     if (in_place)
	  p->p.complex.out = p->p.complex.in;
     else
	  p->p.complex.out = bench_malloc(n * sizeof(bench_complex));

     caset(p->p.complex.out, n, czero);
     caset(p->p.complex.in, n, czero);
     return p;
}

void problem_destroy(struct problem *p)
{
     switch (p->kind) {
	 case PROBLEM_COMPLEX:
	      if (p->p.complex.out != p->p.complex.in)
		   bench_free(p->p.complex.out);
	      bench_free(p->p.complex.in);
	      break;
	 case PROBLEM_REAL:
	      /* TODO */
	      ;
     }

     bench_free(p->n);
     bench_free(p);
}

/* predicates for common cases */
int problem_complex_power_of_two(struct problem *p, int in_place)
{
     return (p->kind == PROBLEM_COMPLEX &&
	     p->rank == 1 &&
	     power_of_two(p->n[0]) &&
	     (in_place ? (p->p.complex.in == p->p.complex.out)
	      : (p->p.complex.in != p->p.complex.out)));
}
