#include "bench-user.h"
#include <math.h>

BEGIN_BENCH_DOC
BENCH_DOC("name", "minfft")
BENCH_DOC("version", "1.2.2")
BENCH_DOC("year", "2022")
BENCH_DOC("author", "Alexander Mukhin")
BENCH_DOC("language", "C")
BENCH_DOC("email", "alexander.i.mukhin@gmail.com")
BENCH_DOC("url", "https://github.com/aimukhin/minfft/")
BENCH_DOC("url-was-valid-on", "Sat Jan  8 16:52:09 UTC 2022")
BENCH_DOC("copyright", "MIT License")
END_BENCH_DOC

#include "minfft.h"

minfft_aux *a; // aux data

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_unpacked(p,out,-1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_unpacked(p,in,-1.0);
}

int can_do(struct problem *p)
{
	unsigned int i;
	for (i=0; i<p->rank; ++i)
		if (!power_of_two(p->n[i]))
			return 0;
	return 1;
}

void setup(struct problem *p)
{
	BENCH_ASSERT(can_do(p));
	if (p->kind==PROBLEM_COMPLEX)
		a=minfft_mkaux_dft(p->rank,(int*)p->n);
	else
		a=minfft_mkaux_realdft(p->rank,(int*)p->n);
	BENCH_ASSERT(a!=NULL);
}

void doit(int iter, struct problem *p)
{
	int i;
	if (p->kind==PROBLEM_COMPLEX)
		// complex DFT
		if (p->sign<0)
			// forward transform
			for (i=0; i<iter; ++i)
				minfft_dft(p->in,p->out,a);
		else
			// inverse transform
			for (i=0; i<iter; ++i)
				minfft_invdft(p->in,p->out,a);
	else
		// real DFT
		if (p->sign<0)
			// forward transform
			for (i=0; i<iter; ++i)
				minfft_realdft(p->in,p->out,a);
		else
			// inverse transform
			for (i=0; i<iter; ++i)
				minfft_invrealdft(p->in,p->out,a);
}

void done(struct problem *p)
{
	UNUSED(p);
	minfft_free_aux(a);
}
