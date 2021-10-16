/* this program is in the public domain */

#include <iostream>
#include "bench-user.h"
#include "ducc0/infra/threading.cc"
#include "ducc0/fft/fft.h"
#include <cmath>
#include <cstring>

using namespace std;
using namespace ducc0;

BEGIN_BENCH_DOC
BENCH_DOC("name", "duccfft")
BENCH_DOC("author", "Martin Reinecke")
BENCH_DOC("year", "2021")
BENCH_DOC("version", "1.0")
BENCH_DOC("language", "C++")
BENCH_DOC("url", "https://gitlab.mpcdf.mpg.de/mtr/ducc")
BENCH_DOC("url-was-valid-on", "Fri Jul 23 23:06:24 ACST 2020")
BENCH_DOC("copyright", "GPLv2+")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return true;
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_1d_fftpack(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_1d_fftpack(p, in, -1.0);
}

static fmav_info::shape_t axes;
static cfmav<complex<bench_real>> in_c;
static vfmav<complex<bench_real>> out_c;
static cfmav<bench_real> in_r;
static vfmav<bench_real> out_r;

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));

     fmav_info::shape_t shape(p->rank);
     axes.resize(p->rank);
     for (int i=0; i<p->rank; ++i) {
        shape[i] = p->n[i];
        axes[i] = i;
        }

     if (p->kind == PROBLEM_COMPLEX) {
        auto in = reinterpret_cast<complex<bench_real> *>(p->in);
        auto out = reinterpret_cast<complex<bench_real> *>(p->out);
        in_c.assign(cfmav<complex<bench_real>>(in, shape));
        out_c.assign(vfmav<complex<bench_real>>(out, shape));
     } else {
        auto in = reinterpret_cast<bench_real *>(p->in);
        auto out = reinterpret_cast<bench_real *>(p->out);
        in_r.assign(cfmav<bench_real>(in, shape));
        out_r.assign(vfmav<bench_real>(out, shape));
     }
}

void doit(int iter, struct problem *p)
{
     if (p->kind == PROBLEM_COMPLEX) {
        for (int i = 0; i < iter; ++i) {
          c2c(in_c,out_c,axes,p->sign==-1,bench_real(1));
	       }
     } else {
	       for (int i = 0; i < iter; ++i) {
          r2r_fftpack(in_r,out_r,axes,p->sign==-1,p->sign==-1,bench_real(1));
	       }
	    }
}

void done(struct problem *p)
{
}
