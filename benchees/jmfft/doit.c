/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "jmfftc")
BENCH_DOC("package", "JMFFT")
BENCH_DOC("author", "Jean-Marie Teuler")
BENCH_DOC("year", "2000")
BENCH_DOC("version", "1.0")
BENCH_DOC("language", "C")
BENCH_DOC("email", "teuler@virgo.infn.it")
BENCH_DOC("email", "Jalel.Chergui@idris.fr")
BENCH_DOC("url", "http://www.idris.fr/data/publications/JMFFT/")
BENCH_DOC("url-was-valid-on", "Sat Aug 31 02:45:22 EDT 2002")
BENCH_DOC("notes", "Designed to emulate the FFTs in Cray's SCILIB.")
BENCH_DOC("notes", "A Fortran-90 version is also available.")
BENCH_DOC("copyright",
"The JMFFT web page says that the library is distributed under the terms "
"of the GNU General Public License as published by the Free Software "
"Foundation (version 2 or later).  However, the source package itself "
"contains no mention of the GPL, or copyright notices per se.  Several "
"of the source files, on the other hand, contain the following notice:\n"
"\n"
"JMFFTLIB : A library of portable fourier transform subroutines "
"emulating Cray SciLib\n"
"Author   : Jean-Marie Teuler, CNRS-IDRIS (teuler@idris.fr)\n"
"\n"
"Permission is granted to copy and distribute this file or modified "
"versions of this file for no fee, provided the copyright notice and "
"this permission notice are preserved on all copies.")
END_BENCH_DOC

#include "jmfft.h"

REAL8 *table = 0, *work = 0;

int can_do(struct problem *p)
{
     return (DOUBLE_PRECISION &&
	     p->rank >= 1 && p->rank <= 3 &&
	     (p->rank != 1 || p->n[0] >= (p->kind == PROBLEM_REAL ? 3 : 2))
	     && (p->kind == PROBLEM_COMPLEX || p->n[p->rank - 1] % 2 == 0));
}

#define MAX2(a,b) ((a) > (b) ? (a) : (b))
#define MAX3(a,b,c) MAX2(a, MAX2(b,c))
#define D2(n) ((n)/2 + 1)
#define DD2(ip, n) ((ip) ? 2 * D2(n) : (n))

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_unpacked(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_unpacked(p, in, -1.0);
}

void copy_r2c(struct problem *p, bench_complex *out)
{
     if (problem_in_place(p) || p->sign > 0)
          copy_r2c_unpacked(p, out);
     else
          copy_r2c_packed(p, out);
}

void copy_c2r(struct problem *p, bench_complex *in)
{
     if (problem_in_place(p) || p->sign > 0)
          copy_c2r_unpacked(p, in);
     else
          copy_c2r_packed(p, in);
}

void setup(struct problem *p)
{
     REAL8 *x = (REAL8 *) p->in;
     REAL8 *y = (REAL8 *) p->out;
     int ip = problem_in_place(p);
     if (p->kind == PROBLEM_COMPLEX) {
	  switch (p->rank) {
	      case 1:
		   table = (REAL8*)bench_malloc(sizeof(REAL8)*(100+8*p->n[0]));
		   work = (REAL8*)bench_malloc(sizeof(REAL8) * 8*p->n[0]);
		   ccfft(0, p->n[0], 1.0, x, y, table, work, 0);
		   break;
	      case 2:
		   table = (REAL8*)bench_malloc(sizeof(REAL8) * 
						(100 + 2*(p->n[0]+p->n[1])));
		   work = (REAL8*)bench_malloc(sizeof(REAL8) * 512
					       * MAX2(p->n[0], p->n[1]));
		   ccfft2d(0, p->n[1], p->n[0], 1.0, x, p->n[1], y, p->n[1],
			   table, work, 0);
		   break;
	      case 3:
		   table = (REAL8*)bench_malloc(sizeof(REAL8) * 
						(100 + 
						 2*(p->n[0]+p->n[1]+p->n[2])));
		   work = (REAL8*)bench_malloc(sizeof(REAL8) * 512
					       * MAX3(p->n[0], p->n[1],
						      p->n[2]));
		   ccfft3d(0, p->n[2], p->n[1], p->n[0], 1.0, 
			   x, p->n[2], p->n[1], y, p->n[2], p->n[1],
			   table, work, 0);
		   break;
	      default:
		   BENCH_ASSERT(0);
	  }
     }
     else {
	  switch (p->rank) {
	      case 1:
		   table = (REAL8*)bench_malloc(sizeof(REAL8)*(100+4*p->n[0]));
		   work = (REAL8*)bench_malloc(sizeof(REAL8) * (4+4*p->n[0]));
		   if (p->sign < 0)
			scfft(0, p->n[0], 1.0, x, y, table, work, 0);
		   else
			csfft(0, p->n[0], 1.0, x, y, table, work, 0);
		   break;
	      case 2:
		   table = (REAL8*)bench_malloc(sizeof(REAL8) * 
						(100 + 2*(p->n[0]+p->n[1])));
		   work = (REAL8*)bench_malloc(sizeof(REAL8) * 512
					       * MAX2(p->n[0], p->n[1]));
		   
		   if (p->sign < 0)
			scfft2d(0, p->n[1], p->n[0], 1.0,
				x, DD2(ip, p->n[1]), y, D2(p->n[1]),
				table, work, 0);
		   else
			csfft2d(0, p->n[1], p->n[0], 1.0,
				x, D2(p->n[1]), y, DD2(1, p->n[1]),
				table, work, 0);
		   break;
	      case 3:
		   table = (REAL8*)bench_malloc(sizeof(REAL8) * 
						(100 + 
						 2*(p->n[0]+p->n[1]+p->n[2])));
		   work = (REAL8*)bench_malloc(sizeof(REAL8) * 512
					       * MAX3(p->n[0], p->n[1],
						      p->n[2]));
		   if (p->sign < 0)
			scfft3d(0, p->n[2], p->n[1], p->n[0], 1.0,
				x, DD2(ip, p->n[2]), p->n[1],
				y, D2(p->n[2]), p->n[1],
				table, work, 0);
		   else
			csfft3d(0, p->n[2], p->n[1], p->n[0], 1.0,
				x, D2(p->n[2]), p->n[1],
				y, DD2(1, p->n[2]), p->n[1],
				table, work, 0);
		   break;
	      default:
		   BENCH_ASSERT(0);
	  }
     }
}

void doit(int iter, struct problem *p)
{
     int n0 = p->n[0];
     int isign = p->sign;
     int ip = problem_in_place(p);
     REAL8 *x = (REAL8 *) p->in;
     REAL8 *y = (REAL8 *) p->out;
     int i;

     if (p->kind == PROBLEM_COMPLEX) {
	  switch (p->rank) {
	      case 1:
		   for (i = 0; i < iter; ++i)
			ccfft(isign, n0, 1.0, x, y, table, work, 0);
		   break;
	      case 2: {
		   int n1 = p->n[1];
		   for (i = 0; i < iter; ++i)
			ccfft2d(isign, n1,n0, 1.0, x,n1, y,n1, table, work, 0);
		   break;
	      }
	      case 3: {
		   int n2 = p->n[2], n1 = p->n[1];
		   for (i = 0; i < iter; ++i)
			ccfft3d(isign, n2, n1, n0, 1.0, 
				x,n2,n1, y,n2,n1, table, work, 0);
		   break;
	      }
	      default:
		   BENCH_ASSERT(0);
	  }
     }
     else {
	  switch (p->rank) {
	      case 1:
		   if (isign < 0)
			for (i = 0; i < iter; ++i)
			     scfft(isign, n0, 1.0, x, y, table, work, 0);
		   else
			for (i = 0; i < iter; ++i)
			     csfft(isign, n0, 1.0, x, y, table, work, 0);
		   break;
	      case 2: {
		   int n1 = p->n[1];
		   int ld = DD2(ip || isign > 0, n1), ld2 = D2(n1);
		   if (isign < 0)
			for (i = 0; i < iter; ++i)
			     scfft2d(isign,n1,n0,1.0,x,ld,y,ld2, table,work,0);
		   else
			for (i = 0; i < iter; ++i)
			     csfft2d(isign,n1,n0,1.0,x,ld2,y,ld, table,work,0);
		   break;
	      }
	      case 3: {
		   int n2 = p->n[2], n1 = p->n[1];
		   int ld = DD2(ip || isign > 0, n2), ld2 = D2(n2);
		   if (isign < 0)
			for (i = 0; i < iter; ++i)
			     scfft3d(isign, n2,n1,n0, 1.0, 
				     x,ld,n1, y,ld2,n1, table,work,0);
		   else
			for (i = 0; i < iter; ++i)
			     csfft3d(isign, n2,n1,n0, 1.0, 
				     x,ld2,n1, y,ld,n1, table,work,0);
		   break;
	      }
	      default:
		   BENCH_ASSERT(0);
	  }
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     bench_free(table);
     bench_free(work);
}
