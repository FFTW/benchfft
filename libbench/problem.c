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

/* $Id: problem.c,v 1.6 2001-07-07 21:56:10 athena Exp $ */

#include "config.h"
#include "bench.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>


/* parse a problem description, return a problem */
/* TODO: only knows about complex for now */
struct problem *problem_parse(const char *s)
{
     int n;
     int in_place = 0;
     struct problem *p;

     p = bench_malloc(sizeof(struct problem));

     p->kind = PROBLEM_COMPLEX;
     p->sign = -1;
     p->in = 0;
     p->out = 0;
     p->userinfo = 0;

     for (;;) {
	  switch (tolower(*s)) {
	      case 'i': in_place = 1; ++s; continue;
	      case 'o': in_place = 0; ++s; continue;
	      case 'f': 
	      case '-': p->sign = -1; ++s; continue;
	      case 'b': 
	      case '+': p->sign = 1; ++s; continue;
	      case 'r': p->kind = PROBLEM_REAL; ++s; continue;
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

     problem_alloc(p, in_place);
     return p;
}

void problem_destroy(struct problem *p)
{
     BENCH_ASSERT(p);
     problem_free(p);
     bench_free(p);
}

/* predicates for common cases */
int problem_power_of_two(struct problem *p, int in_place)
{
     unsigned int i;

     for (i = 0; i < p->rank; ++i)
	  if (!(power_of_two(p->n[i])))
	       return 0;

     return (in_place ? (p->in == p->out) : (p->in != p->out));
}

int problem_complex_power_of_two(struct problem *p, int in_place)
{
     if (p->kind != PROBLEM_COMPLEX)
	  return 0;

     return problem_power_of_two(p, in_place);
}

int problem_real_power_of_two(struct problem *p, int in_place)
{
     if (p->kind != PROBLEM_REAL)
	  return 0;

     return problem_power_of_two(p, in_place);
}
