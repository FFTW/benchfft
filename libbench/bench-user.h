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

/* $Id: bench-user.h,v 1.2 2001-07-05 14:27:19 athena Exp $ */
#ifndef __BENCH_USER_H__
#define __BENCH_USER_H__

/* benchmark program definitions for user code */
#include "config.h"

#if HAVE_STDDEF_H
#include <stddef.h>
#endif

typedef double bench_real;

typedef struct {
     bench_real re, im;
} bench_complex;

#define c_re(c)  ((c).re)
#define c_im(c)  ((c).im)

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

#ifdef CC
#define CC_DOC BENCH_DOC("cc", CC)
#else
#define CC_DOC /* none */
#endif

#ifdef CFLAGS
#define CFLAGS_DOC BENCH_DOC("cflags", CFLAGS)
#else
#define CFLAGS_DOC /* none */
#endif

#ifdef F77
#define F77_DOC BENCH_DOC("f77", F77)
#else
#define F77_DOC /* none */
#endif

#ifdef FFLAGS
#define FFLAGS_DOC BENCH_DOC("fflags", FFLAGS)
#else
#define FFLAGS_DOC /* none */
#endif

#define BEGIN_BENCH_DOC				\
const struct bench_doc bench_doc[] = {		\
    CC_DOC					\
    CFLAGS_DOC					\
    F77_DOC					\
    FFLAGS_DOC

#define BENCH_DOC(key, val) { key, val },

#define END_BENCH_DOC				\
     {0, 0}};
    
#endif /* __BENCH_USER_H__ */
