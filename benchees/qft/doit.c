/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "qft")
BENCH_DOC("author", "Gary A. Sitton")
BENCH_DOC("year", "1994")
BENCH_DOC("language", "C")
BENCH_DOC("email", "gary.sitton@wg.waii.com")
BENCH_DOC("bibitem", 
	  "H. Guo, G. A. Sitton, C. S. Burrus, The Quick Discrete Fourier Transform, Proc. ICASSP, April 1994.")
BENCH_DOC("copyright", "As I developed this software in conjunction with "
	  "Sid Burrus and Haitao Guo at Rice University, I "
	  "freely put this software in the public domain. "
	  "As such, it requires no banner preservation, "
	  "restrictions, or acknowledgements other than "
	  "those which the user wishes to add.")
BENCH_DOC("notes", "Received in personal communication with the author.")
BENCH_DOC("notes", "Complex data are stored in separate real/imag arrays.")
BENCH_DOC("notes", "inverse real transform is scaled")
END_BENCH_DOC

#include "qftsubs.h"

int can_do(struct problem *p)
{
     return (SINGLE_PRECISION &&
	     p->rank == 1 &&
	     problem_power_of_two(p, 0) &&
	     p->n[0] >= NN);
}

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
     BENCH_ASSERT(qftinit(p->n[0]));
}

void doit(int iter, struct problem *p)
{
     int i, n = p->n[0];

     if (p->kind == PROBLEM_COMPLEX) {
	  float *rin = (float*) p->in, *iin;
	  float *rout = (float*) p->out, *iout;
	  iin = rin + n; iout = rout + n;

	  if (p->sign < 0)
	       for (i = 0; i < iter; ++i)
		    cfqft(rout, iout, rin, iin, n, n);
	  else
	       for (i = 0; i < iter; ++i)
		    cfqft(iout, rout, iin, rin, n, n);
     }
     else if (p->kind == PROBLEM_REAL) {
	  float *rin = (float*) p->in, *iin;
	  float *rout = (float*) p->out, *iout;
	  iin = rin + n/2 + 1; iout = rout + n/2 + 1;

	  if (p->sign < 0)
	       for (i = 0; i < iter; ++i)
		    rfqft(rout, iout, rin, n, n);
	  else
	       for (i = 0; i < iter; ++i)
		    riqft(rin, iin, rout, n, n);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     qftclos();
}

void copy_c2c_from(struct problem *p, bench_complex *in)
{
     bench_real *x = (bench_real *) p->in;
     copy_c2ri(in, x, x + p->n[0], p->n[0]);
}

void copy_c2c_to(struct problem *p, bench_complex *out)
{
     bench_real *x = (bench_real *) p->out;
     copy_ri2c(x, x + p->n[0], out, p->n[0]);
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_1d_unpacked_ri(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_unpacked_ri(p, in, -1.0);
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     if (p->kind == PROBLEM_REAL)
	  unnormalize(p, out, 1);
}
