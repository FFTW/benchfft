if HAVE_F77
noinst_PROGRAMS=doit
noinst_LIBRARIES=libfft.a
endif

INCLUDES = -I${top_srcdir}/libbench

doit_SOURCES=doit.c
doit_LDADD=libfft.a ${top_builddir}/libbench/libbench.a

# use separate library so that automake links doit as a C program and
# not FORTRAN
libfft_a_SOURCES=fft.f