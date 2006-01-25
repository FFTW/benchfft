/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "cwplib")
BENCH_DOC("author", "Dave Hale")
BENCH_DOC("year", "1998")
BENCH_DOC("version", "Release 39")
BENCH_DOC("language", "C")
BENCH_DOC("url", "http://www.cwp.mines.edu/cwpcodes/")
BENCH_DOC("url-was-valid-on", "Tue Jan 24 18:54:17 EST 2006")
BENCH_DOC("copyright", "Copyright (c) 2005, Colorado School of Mines, all rights reserved.\n"
"\n"
"Warranty Disclaimer:\n"
"NO GUARANTEES, OR WARRANTIES, EITHER EXPRESS OR IMPLIED, ARE PROVIDED BY \n"
"CWP, CSM, ANY EMPLOYEE OR MEMBER OF THE AFORESAID ORGANIZATIONS, OR BY \n"
"ANY CONTRIBUTOR TO THIS SOFTWARE PACKAGE, REGARDING THE ACCURACY, SAFETY, \n"
"FITNESS FOR A PARTICULAR PURPOSE, OR ANY OTHER QUALITY OR CHARACTERISTIC \n"
"OF THIS SOFTWARE, OR ANY DERIVATIVE SOFTWARE.\n"
"\n"
"Export Restriction Disclaimer:\n"
"We believe that CWP/SU: Seismic Un*x is a low technology product that does\n"
"not appear on the Department of Commerce CCL list of restricted exports.\n"
"Accordingly, we believe that our product meets the qualifications of\n"
"an ECCN (export control classification number) of EAR99 and we believe\n"
"it fits the qualifications of NRR (no restrictions required), and\n"
"is thus not subject to export restrictions of any variety.\n"
"\n"
"Limited License:\n"
"The CWP/SU Seismic Un*x package (SU) is not public domain software, but \n"
"it is available free under the following terms and conditions:\n"
"\n"
"1. Permission to use, copy, and modify this software for any purpose without\n"
"fee and within the guidelines set forth below is hereby granted, provided \n"
"that the above copyright notice, the warranty disclaimer, and this \n"
"permission notice appear in all copies, and the name of the \n"
"Colorado School of Mines (CSM) not be used\n"
"in advertising or publicity pertaining to this software without the \n"
"specific, written permission of CSM.\n"
"\n"
"2. The simple repackaging and selling of the SU package as is, as a \n"
"commercial software product, is expressly forbidden without the prior \n"
"written permission of CSM.  Any approved repackaging arrangement will \n"
"carry the following restriction: only a modest profit over reproduction \n"
"costs may be realized by the reproducer.\n")
BENCH_DOC("notes", "Revised by Baoniu Han to handle double precision. 12/14/98")
BENCH_DOC("notes", "Benchmark uses the files pfafft.c, dpfafft.c, and cwp.h, which are included as a part of the CWP/SU Seismic Un*x package.")
END_BENCH_DOC

#include "cwp.h"

/* cwp double-precision routines have _d suffix */
#ifdef BENCHFFT_SINGLE
#  define F(x) x
#else
#  define F(x) x ## _d
#endif

int can_do(struct problem *p)
{
     unsigned int i;

     if (p->kind == PROBLEM_COMPLEX) {
	  /* 
	   * CWPLIB does not really provide a 2D transform.
	   * Nevertheless, source code comments provide instructions
	   * on how to use CWPLIB to compute a 2D transform, which we
	   * take as an implementation.
	   */
	  if (p->rank != 1 && p->rank != 2)
	       return 0;

	  /* must be in place */
	  if (!problem_in_place(p))
	       return 0;

	  for (i = 0; i < p->rank; ++i) {
	       if ((unsigned int) F(npfa)(p->n[i]) != p->n[i])
		    return 0;
	  } 
     } else {
	  /* no 2D instructions are provided */
	  if (p->rank != 1)
	       return 0;

	  /* can be either in place or out of place */

	  if ((unsigned int) F(npfar)(p->n[0]) != p->n[0])
	       return 0;
     }

     return 1;
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_unpacked(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_unpacked(p, in, -1.0);
}

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
     /* nothing to do */
}

void doit(int iter, struct problem *p)
{
     int i;
     int sign = p->sign;
     void *in = p->in;

     if (p->kind == PROBLEM_COMPLEX) {
	  switch (p->rank) {
	      case 1:
	      {
		   int n = p->n[0];
		   for (i = 0; i < iter; ++i) 
			F(pfacc)(sign, n, in);
		   break;
	      }
	      case 2:
	      {
		   int n0 = p->n[0];
		   int n1 = p->n[1];
		   for (i = 0; i < iter; ++i) { 
			F(pfa2cc)(sign, 1, n1, n0, in);
			F(pfa2cc)(sign, 2, n1, n0, in);
		   }
		   break;
	      }
     
	      default:
		   BENCH_ASSERT(0);
	  }
     } else {
	  /* p->rank == 1*/
	  int n = p->n[0];
	  void *out = p->out;

	  if (sign == -1) {
	       /* real->complex */
	       for (i = 0; i < iter; ++i) 
		    F(pfarc)(sign, n, in, out);
	  } else {
	       for (i = 0; i < iter; ++i) 
		    F(pfacr)(sign, n, in, out);
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
}
