/* this program is in the public domain */

#include "bench-user.h"
#include "fourier.h"

BEGIN_BENCH_DOC
BENCH_DOC("name", "cross")
BENCH_DOC("author", "Don Cross")
BENCH_DOC("year", "1998")
BENCH_DOC("language", "C")
BENCH_DOC("email", "dcross@intersrv.com")
BENCH_DOC("url", "http://www.intersrv.com/~dcross/fft.html")
BENCH_DOC("url-was-valid-on", "Mon Sep  2 21:54:04 EDT 2002")
BENCH_DOC("notes", "Complex data are stored in separate real/imag arrays.")
BENCH_DOC("notes", "The forward transform is scaled.")
BENCH_DOC("copyright",
"                         *** SMALL REQUESTS ****\n"
"\n"
"If you want to give away copies of this source code, that's fine, so long\n"
"as you do the following:\n"
"\n"
"- Do not charge any money for this source code, except for possibly a\n"
"  reasonable fee to cover postage, disk duplication, etc.  I wrote this\n"
"  code and I want it to be FREE to EVERYONE!\n"
"\n"
"- Do not remove my name, e-mail address, or URL from any of the files in\n"
"  this collection.\n"
"\n"
"- Please keep this readme.txt file with the source and headers so that others\n"
"  can get in touch with me to ask questions and/or find my web page to read\n"
"  the online tutorial.\n"
"\n"
"- If you make any modifications to the source code, please add comments to\n"
"  it which include your name, e-mail address, web page URL (if any), and\n"
"  explain what you did to the code.\n"
"\n"
"- If you use this source code in an interesting program, please let me know.\n"
"  I promise will never try to get money from you, even if you use this code\n"
"  in a commercial program.  I just want to know what kind of clever and\n"
"  creative things people do with Fourier Transforms.\n"
)
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (p->rank == 1 &&
	     p->kind == PROBLEM_COMPLEX &&
	     problem_power_of_two(p, 0));
}

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
}

void copy_c2c_from(struct problem *p, bench_complex *in)
{
     bench_real *x = (bench_real *) p->in;
     copy_c2ri(in, x, x + p->size, p->size);
}

void copy_c2c_to(struct problem *p, bench_complex *out)
{
     bench_real *x = (bench_real *) p->out;
     copy_ri2c(x, x + p->size, out, p->size);
}

void after_problem_ccopy_to(struct problem *p, bench_complex *out)
{
     unnormalize(p, out, -1);
}

void doit(int iter, struct problem *p)
{
     int n = p->n[0], dir = p->sign > 0 ? 0 : 1;
     bench_real *ri = (bench_real *) p->in, *ii;
     bench_real *ro = (bench_real *) p->out, *io;
     int i;
     ii = ri + n; io = ro + n;

#ifdef BENCHFFT_SINGLE
     for (i = 0; i < iter; ++i)
	  fft_float(n, dir, ri, ii, ro, io);
#else
     for (i = 0; i < iter; ++i)
	  fft_double(n, dir, ri, ii, ro, io);
#endif
}

void done(struct problem *p)
{
     UNUSED(p);
}
