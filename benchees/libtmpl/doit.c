/* this program is in the public domain */

#include "bench-user.h"
#include <libtmpl/include/tmpl_fft.h>
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "libtmpl")
BENCH_DOC("author", "Ryan Maguire")
BENCH_DOC("year", "2023")
BENCH_DOC("version", "(none)")
BENCH_DOC("language", "C")
BENCH_DOC("url", "https://github.com/ryanmaguire/libtmpl")
BENCH_DOC("url-was-valid-on", "Tue Mar 18 14:09:12 EDT 2025")
BENCH_DOC("copyright", "GPL version 3 or later")
END_BENCH_DOC

int can_do(struct problem *p)
{
	return (DOUBLE_PRECISION && p->rank == 1 && !problem_in_place(p) && p->kind == PROBLEM_COMPLEX);
}

void setup(struct problem *p)
{
	BENCH_ASSERT(can_do(p));
}

void doit(int iter, struct problem *p)
{
	int i;
	size_t n = (size_t) p->n[0];
	tmpl_ComplexDouble *in = (tmpl_ComplexDouble*) p->in, *out = (tmpl_ComplexDouble*) p->out;

	if (p->sign == -1) {
		for (i = 0; i < iter; ++i) {
			tmpl_CDouble_FFT(in, out, n);
		}
	} else {
		for (i = 0; i < iter; ++i) {
			tmpl_CDouble_IFFT(in, out, n);
		}
	}
}

void done(struct problem *p)
{
	BENCH_ASSERT(can_do(p));
}
