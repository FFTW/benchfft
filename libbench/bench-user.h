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

/* $Id: bench-user.h,v 1.1 2001-07-04 22:50:36 athena Exp $ */

/* benchmark program definitions for user code */
#include "config.h"

#if HAVE_STDDEF_H
#include <stddef.h>
#endif

typedef double bench_real;

typedef struct {
     bench_real re, im;
} bench_complex;

struct problem {
     enum { PROBLEM_COMPLEX, PROBLEM_REAL } kind;
     int rank;
     int *n;  
     int size;  /* total size of input = PROD n[i] */

     void *userinfo; /* user can store whatever */
     union {
	  struct {
	       int sign;
	       bench_complex *in;
	       bench_complex *out;
	  } complex;

	  struct {
	       /* TODO */
	       int dummy; /* prevent compiler warning */
	  } real;
     } p;
};

extern int can_do(struct problem *p);
extern void setup(struct problem *p);
extern void doit(int iter, struct problem *p);
extern void done(struct problem *p);

#define power_of_two(n) (((n) > 0) && (((n) & ((n) - 1)) == 0))
int problem_complex_power_of_two(struct problem *p, int in_place);

/**************************************************************
 * malloc
 **************************************************************/
extern void *bench_malloc(size_t size);
extern void bench_free(void *ptr);

/**************************************************************
 * alloca
 **************************************************************/
#ifdef HAVE_ALLOCA_H
#include <alloca.h>
#endif

#ifdef HAVE_ALLOCA
/* use alloca if available */
#define STACK_MALLOC(x) alloca(x)
#define STACK_FREE(x) 

#else
/* use malloc instead of alloca */
#define STACK_MALLOC(x) bench_malloc(x)
#define STACK_FREE(x) bench_free(x)
#endif


/**************************************************************
 * assert
 **************************************************************/
extern void bench_assertion_failed(const char *s, int line, char *file);
#define BENCH_ASSERT(ex)						 \
      (void)((ex) || (bench_assertion_failed(#ex, __LINE__, __FILE__), 0))

#define UNUSED(x) (void)x

/***************************************
 * Documentation strings
 ***************************************/
struct bench_doc {
     const char *key;
     const char *val;
};

extern const struct bench_doc bench_doc[];

#define BEGIN_BENCH_DOC				\
const struct bench_doc bench_doc[] = {

#define BENCH_DOC(key, val) { key, val },

#define END_BENCH_DOC				\
     {0, 0}};
    
