/* this program is in the public domain */

#include "bench-user.h"
#include "pocketfft_hdronly.h"
#include <cmath>
#include <cstring>

using namespace std;
using namespace pocketfft;

BEGIN_BENCH_DOC
BENCH_DOC("name", "pocketfft_cxx")
BENCH_DOC("author", "Martin Reinecke")
BENCH_DOC("year", "2019")
BENCH_DOC("version", "1.0")
BENCH_DOC("language", "C++")
BENCH_DOC("url", "https://gitlab.mpcdf.mpg.de/mtr/pypocketfft")
BENCH_DOC("url-was-valid-on", "Fri Jul 23 23:06:24 ACST 2020")
BENCH_DOC("copyright", "3 clause BSDL")
END_BENCH_DOC

int can_do(struct problem * /*p*/)
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

static shape_t shape, axes;
static stride_t strides;

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
     
     shape.resize(p->rank);
     strides.resize(p->rank);
     axes.resize(p->rank);
     for (int i=0; i<p->rank; ++i) {
       shape[i] = p->n[i];
       axes[i] = i;
     }
     ptrdiff_t str=sizeof(bench_real);
     if (p->kind == PROBLEM_COMPLEX) str*=2;
     for (int i=p->rank-1; i>=0; --i) {
       strides[i] = str;
       str *= shape[i];
     }
}

void doit(int iter, struct problem *p)
{
     if (p->kind == PROBLEM_COMPLEX) {
        auto in = reinterpret_cast<complex<bench_real> *>(p->in);
        auto out = reinterpret_cast<complex<bench_real> *>(p->out);
        for (int i = 0; i < iter; ++i) {
          c2c(shape, strides, strides, axes,p->sign==-1,in, out,bench_real(1));
	       }
     } else {
        auto in = reinterpret_cast<bench_real *>(p->in);
        auto out = reinterpret_cast<bench_real *>(p->out);
	       for (int i = 0; i < iter; ++i) {
          r2r_fftpack(shape, strides,strides,axes,p->sign==-1,p->sign==-1,in,out,bench_real(1));
	       }
	    }
}

void done(struct problem *p)
{
}
