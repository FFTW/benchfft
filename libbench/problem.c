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

/* $Id: problem.c,v 1.2 2001-07-05 16:49:43 athena Exp $ */

#include "config.h"
#include "bench.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>


/* parse a problem description, return a problem */
/* TODO: only knows complex 1D for now */
struct problem *problem_parse(const char *s)
{
     int n;
     int in_place = 0;
     int sign = -1;
     struct problem *p;
     const bench_complex czero = {0, 0};

     p = bench_malloc(sizeof(struct problem));

     p->kind = PROBLEM_COMPLEX;
     p->userinfo = 0;

     for (;;) {
	  switch (*s) {
	      case 'i':
	      case 'I':
		   in_place = 1;
		   ++s;
		   continue;
	      case 'f':
	      case 'F':
		   sign = -1;
		   ++s;
		   continue;
	      case 'b':
	      case 'B':
		   sign = 1;
		   ++s;
		   continue;
	      default:
		   ;
	  }
	  break;
     }

     p->rank = 0;
     p->size = 1;      /* the product of 0 things is 1 */

 accept_digit:
     n = 0;

     BENCH_ASSERT(isdigit(*s));

     while (isdigit(*s)) {
	  n = n * 10 + (*s - '0');
	  ++s;
     }

     p->n[p->rank] = n;
     p->size *= n;
     ++p->rank;

     BENCH_ASSERT(p->rank < MAX_RANK);

     if (*s == 'x' || *s == 'X' || *s == '*' || *s == ',') {
	  ++s;
	  goto accept_digit;
     }

     p->p.complex.sign = sign;
     p->p.complex.in = bench_malloc(p->size * sizeof(bench_complex));

     if (in_place)
	  p->p.complex.out = p->p.complex.in;
     else
	  p->p.complex.out = bench_malloc(p->size * sizeof(bench_complex));

     caset(p->p.complex.out, p->size, czero);
     caset(p->p.complex.in, p->size, czero);
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
