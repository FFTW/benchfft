if DOUBLE
if HAVE_F90
PRG=doit
noinst_LIBRARIES=libfft.a
endif
endif

INCLUDES = $(INCLBENCH)
EXTRA_DIST = sgfft.tar.gz

libfft_a_SOURCES=basis_defs.f90 sg_ctrig.f90 sg_fft.f90 sg_fftpad.f90	\
sg_fftpx.f90 sg_fftrisc.f90 sg_fftx.f90 sg_ffty.f90 sg_fftz.f90

doit_SOURCES=doit.c
doit_LDADD=libfft.a $(LIBBENCH) @FLIBS@

include ../Makefile.common
