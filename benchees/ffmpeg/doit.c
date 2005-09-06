/* this program is in the public domain */

#include "bench-user.h"
#include <math.h>

#include "dsputil.h"

BEGIN_BENCH_DOC
BENCH_DOC("name", "ffmpeg")
BENCH_DOC("author", "Fabrice Bellard")
BENCH_DOC("year", "2002")
BENCH_DOC("version", LIBAVCODEC_IDENT)
BENCH_DOC("language", "C")
BENCH_DOC("email", "fabrice@bellard.org")
BENCH_DOC("url", "http://ffmpeg.sourceforge.net/")
BENCH_DOC("url-was-valid-on", "Sun Sep  4 21:02:50 EDT 2005")
BENCH_DOC("notes", "Complex-data FFT routines from libavcodec, with whatever compilation flags were used for the libavcodec installed on the system.")
BENCH_DOC("notes", "Benchmarked ff_fft_permute + ff_fft_calc pair of routines, so that we are benchmarking the normal-order FFT.")
BENCH_DOC("copyright",
"Copyright (c) 2002 Fabrice Bellard.\n"
"\n"
"This library is free software; you can redistribute it and/or\n"
"modify it under the terms of the GNU Lesser General Public\n"
"License as published by the Free Software Foundation; either\n"
"version 2 of the License, or (at your option) any later version.\n"
"\n"
"This library is distributed in the hope that it will be useful,\n"
"but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU\n"
"Lesser General Public License for more details.\n"
"\n"
"You should have received a copy of the GNU Lesser General Public\n"
"License along with this library; if not, write to the Free Software\n"
"Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA\n")
END_BENCH_DOC

int can_do(struct problem *p)
{
     return (SINGLE_PRECISION &&
	     p->rank == 1 &&
	     p->kind == PROBLEM_COMPLEX &&
	     p->n[0] >= 8 && /* 20050313 version crashes for smaller sizes */
	     p->n[0] <= 65536 && /* 20050313 has a uint16_t bitrev table */
	     problem_power_of_two(p, 1));
}

FFTContext cntxt;

void setup(struct problem *p)
{
     BENCH_ASSERT(can_do(p));
     ff_fft_init(&cntxt, log_2(p->n[0]), p->sign > 0);
}

void doit(int iter, struct problem *p)
{
     FFTComplex *x = (FFTComplex *) p->in;
     int i;

     for (i = 0; i < iter; ++i) {
	  ff_fft_permute(&cntxt, x);
	  ff_fft_calc(&cntxt, x);
     }
}

void done(struct problem *p)
{
     UNUSED(p);
     ff_fft_end(&cntxt);
}
