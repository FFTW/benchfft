if HAVE_F77
noinst_PROGRAMS=doit
noinst_LIBRARIES=libfft.a
endif

INCLUDES = -I${top_srcdir}/libbench

doit_SOURCES=doit.c
doit_LDADD=libfft.a ${top_builddir}/libbench/libbench.a @FLIBS@

libfft_a_SOURCES=fft.f
