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

/* $Id: zero.c,v 1.5 2002-08-16 22:23:39 athena Exp $ */

#include "config.h"
#include "bench.h"

/* set I/O arrays to zero.  Default routine */
void problem_zero(struct problem *p)
{
     if (p->kind == PROBLEM_COMPLEX) {
	  const bench_complex czero = {0, 0};
	  caset(p->out, p->phys_size * p->batch, czero);
          if (p->in != p->out)
              caset(p->in, p->phys_size * p->batch, czero);
     } else {
	  aset(p->out, p->phys_size * p->batch, 0.0);
          if (p->in != p->out)
              aset(p->in, p->phys_size * p->batch, 0.0);
     }
}
