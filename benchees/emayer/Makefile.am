if DOUBLE
if HAVE_FC
PRG=doit
endif
endif

EXTRA_DIST=fft_radix8_fwd_bench_rc.f90 fft_radix8_rev_bench_rc.f90
INCLUDES = $(INCLBENCH)

doit_SOURCES = doit.c
doit_LDADD = fft_radix8_fwd_bench_rc.o fft_radix8_rev_bench_rc.o	\
	     $(LIBBENCH) $(FCLIBS)

fft_radix8_fwd_bench_rc.o: fft_radix8_fwd_bench_rc.f90
	$(FC) $(FCFLAGS) -c $(FCFLAGS_f90) ${srcdir}/fft_radix8_fwd_bench_rc.f90

fft_radix8_rev_bench_rc.o: fft_radix8_rev_bench_rc.f90
	$(FC) $(FCFLAGS) -c $(FCFLAGS_f90) ${srcdir}/fft_radix8_rev_bench_rc.f90

include ../Makefile.common
