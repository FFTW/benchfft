/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "green-ffts-2.0")
BENCH_DOC("author", "John Green")
BENCH_DOC("email", "green_jt@vsdec.npt.nuwc.navy.mil")
BENCH_DOC("copyright",
	  "This code is public domain, do anything you want to with it.")
BENCH_DOC("language", "C")
BENCH_DOC("notes", "The backward transform is scaled")
END_BENCH_DOC

#include "fftext.h"
#include "fft2d.h"

int can_do(struct problem *p)
{
     return (sizeof(bench_real) == sizeof(float) &&
	     problem_power_of_two(p, 1) &&
	     ((p->kind == PROBLEM_COMPLEX && p->rank >= 1 && p->rank < 4)
	      ||
	      (p->kind == PROBLEM_REAL && p->rank >= 1 && p->rank < 3)));
}


void problem_ccopy_to(struct problem *p, bench_complex *out)
{
     if (p->kind == PROBLEM_COMPLEX) {
	  cacopy(p->out, out, p->size);
	  if (p->sign == 1) { 	  /* undo the scaling */
	       bench_complex x;
	       c_re(x) = p->size;
	       c_im(x) = 0;

	       cascale(out, p->size, x);
	  }
     } else {
	  BENCH_ASSERT(0);
     }
}

void setup(struct problem *p)
{
     switch (p->rank) {
	 case 1: 
	      fftInit(log_2(p->n[0]));
	      break;
	 case 2: 
	      fft2dInit(log_2(p->n[0]), log_2(p->n[1]));
	      break;
	 case 3:
	      fft3dInit(log_2(p->n[0]), log_2(p->n[1]), log_2(p->n[2]));
	      break;
	 default:
	      BENCH_ASSERT(/* can't happen */ 0);
	      break;
     }
}

void doit(int iter, struct problem *p)
{
     int i;
     void *in = p->in;

     if (p->kind == PROBLEM_COMPLEX) {
	  switch (p->rank) {
	      case 1: 
	      {
		   int m0 = log_2(p->n[0]);

		   if (p->sign == -1) {
			for (i = 0; i < iter; ++i) {
			     ffts(in, m0, 1);
			}
		   } else {
			for (i = 0; i < iter; ++i) {
			     iffts(in, m0, 1);
			}
		   }
		   break;
	      }
	      case 2: 
	      {
		   int m0 = log_2(p->n[0]);
		   int m1 = log_2(p->n[1]);

		   if (p->sign == -1) {
			for (i = 0; i < iter; ++i) {
			     fft2d(in, m0, m1);
			}
		   } else {
			for (i = 0; i < iter; ++i) {
			     ifft2d(in, m0, m1);
			}
		   }
		   break;
	      }
	      case 3: 
	      {
		   int m0 = log_2(p->n[0]);
		   int m1 = log_2(p->n[1]);
		   int m2 = log_2(p->n[2]);

		   if (p->sign == -1) {
			for (i = 0; i < iter; ++i) {
			     fft3d(in, m0, m1, m2);
			}
		   } else {
			for (i = 0; i < iter; ++i) {
			     ifft3d(in, m0, m1, m2);
			}
		   }
		   break;
	      }
	  }
     } else {
	  switch (p->rank) {
	      case 1: 
	      {
		   int m0 = log_2(p->n[0]);

		   if (p->sign == -1) {
			for (i = 0; i < iter; ++i) {
			     rffts(in, m0, 1);
			}
		   } else {
			for (i = 0; i < iter; ++i) {
			     riffts(in, m0, 1);
			}
		   }
		   break;
	      }
	      case 2: 
	      {
		   int m0 = log_2(p->n[0]);
		   int m1 = log_2(p->n[1]);

		   if (p->sign == -1) {
			for (i = 0; i < iter; ++i) {
			     rfft2d(in, m0, m1);
			}
		   } else {
			for (i = 0; i < iter; ++i) {
			     rifft2d(in, m0, m1);
			}
		   }
		   break;
	      }
	  }
     }
}

void done(struct problem *p)
{
     switch (p->rank) {
	 case 1: 
	      fftFree();
	      break;
	 case 2: 
	      fft2dFree();
	      break;
	 case 3:
	      fft3dFree();
	      break;
	 default:
	      BENCH_ASSERT(/* can't happen */ 0);
	      break;
     }
}
