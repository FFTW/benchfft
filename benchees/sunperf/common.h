
static const char *mkvers(void)
{
#ifdef HAVE_SUNPERF_VERSION_
     int a, b, c;
     static char buf[128];

     /* 
	using FORTRAN interface.  The documented C interface cannot
	possibly work because it takes a, b, c by value and it returns
	void
     */
     sunperf_version_(&a, &b, &c);
     sprintf(buf, "%d.%d.%d", a, b, c);
     return buf;
#else
     return "unknown";
#endif
}

BEGIN_BENCH_DOC
BENCH_DOC("name", "sunperf")
BENCH_DOCF("version", mkvers)
BENCH_DOC("package", "Sun Performance Library (SUNPERF)")
BENCH_DOC("notes", "Using new fft interface (CFFTC etc.)")
END_BENCH_DOC


void copy_r2c(struct problem *p, bench_complex *out)
{
     copy_r2c_unpacked(p, out);	  
}

void copy_c2r(struct problem *p, bench_complex *in)
{
     copy_c2r_unpacked(p, in);
}

void copy_h2c(struct problem *p, bench_complex *out)
{
     copy_h2c_unpacked(p, out, -1.0);
}

void copy_c2h(struct problem *p, bench_complex *in)
{
     copy_c2h_unpacked(p, in, -1.0);
}

#define MAX2(a,b) ((a) > (b) ? (a) : (b))
#define MAX3(a,b,c) MAX2(a,MAX2(b,c))
