if DOUBLE
if HAVE_F90
PRG=doit
endif
endif

EXTRA_DIST=fft_radix8_fwd_bench_rc.f90 fft_radix8_rev_bench_rc.f90
INCLUDES = $(INCLBENCH)

doit_SOURCES = doit.c
doit_LDADD = fft_radix8_fwd_bench_rc.o fft_radix8_rev_bench_rc.o	\
	     $(LIBBENCH) $(F90LIBS)

fft_radix8_fwd_bench_rc.o: fft_radix8_fwd_bench_rc.f90
	$(F90) $(F90FLAGS) -c ${srcdir}/fft_radix8_fwd_bench_rc.f90

fft_radix8_rev_bench_rc.o: fft_radix8_rev_bench_rc.f90
	$(F90) $(F90FLAGS) -c ${srcdir}/fft_radix8_rev_bench_rc.f90

include ../Makefile.common
