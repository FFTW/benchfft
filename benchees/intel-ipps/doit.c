#include "bench-user.h"
#include <math.h>
#include <stdio.h>
#include <ipps.h>

static const char *mkvers(void)
{
    static char buf[1024];
    const IppLibraryVersion *lib = ippsGetLibVersion();

    sprintf(buf, "%s %s %d.%d.%d.%d", lib->Name, lib->Version,
	    lib->major, lib->minor, lib->majorBuild, lib->build);
    return buf;
}

BEGIN_BENCH_DOC
BENCH_DOC("name", "intel-ipps")
BENCH_DOC("package", "Intel Integrated Performance Primitives, Signal Processing (IPPS)")
BENCH_DOCF("version", mkvers)
BENCH_DOC("notes", "Using ``Perm'' format for real transforms.")
END_BENCH_DOC

#define CONCAT(prefix, name) prefix ## name

#ifdef BENCHFFT_SINGLE
#define MANGLEC(name) CONCAT(name, _32fc)
#define MANGLE(name) CONCAT(name, _32f)
#else
#define MANGLEC(name) CONCAT(name, _64fc)
#define MANGLE(name) CONCAT(name, _64f)
#endif

void copy_h2c(struct problem *p, bench_complex *out)
{
    copy_h2c_1d_packed(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
    copy_c2h_1d_packed(p, in, -1.0);
}

static void *thing;
static void *buffer;

int can_do(struct problem *p)
{
    return (p->rank == 1 &&
	    (problem_power_of_two(p, 0) || problem_power_of_two(p, 1)) );
}

void setup(struct problem *p)
{
    int bufsz;

    if (p->kind == PROBLEM_COMPLEX) {
	MANGLEC(ippsFFTInitAlloc_C)(&thing, 
				    log_2(p->n[0]), 
				    IPP_FFT_NODIV_BY_ANY,
				    ippAlgHintAccurate);
	MANGLEC(ippsFFTGetBufSize_C)(thing, &bufsz);
    } else {
	MANGLE(ippsFFTInitAlloc_R)(&thing, 
				   log_2(p->n[0]), 
				   IPP_FFT_NODIV_BY_ANY,
				   ippAlgHintAccurate);
	MANGLE(ippsFFTGetBufSize_R)(thing, &bufsz);
    }

    if (verbose > 2)
	printf("bufsize = %d\n", bufsz);
    buffer = bench_malloc(bufsz);
}

void doit(int iter, struct problem *p)
{
    int i;

    if (p->kind == PROBLEM_COMPLEX) {
	if (p->sign == -1)
	    for (i = 0; i < iter; ++i) {
		MANGLEC(ippsFFTFwd_CToC)(p->in, p->out, thing, buffer);
	    }
	else 
	    for (i = 0; i < iter; ++i) {
		MANGLEC(ippsFFTInv_CToC)(p->in, p->out, thing, buffer);
	    }
    } else {
	if (p->sign == -1)
	    for (i = 0; i < iter; ++i) {
		MANGLE(ippsFFTFwd_RToPerm)(p->in, p->out, thing, buffer);
	    }
	else 
	    for (i = 0; i < iter; ++i) {
		MANGLE(ippsFFTInv_PermToR)(p->in, p->out, thing, buffer);
	    }
    }
}

void done(struct problem *p)
{
    if (p->kind == PROBLEM_COMPLEX) {
	MANGLEC(ippsFFTFree_C)(thing);
    } else {
	MANGLE(ippsFFTFree_R)(thing);
    }
    thing = 0;
    bench_free(buffer);
}
