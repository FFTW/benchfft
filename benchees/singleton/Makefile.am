if HAVE_F77
noinst_PROGRAMS=doit
endif

INCLUDES = -I${top_srcdir}/libbench

doit_SOURCES=doit.c fft.f
doit_LDADD=${top_builddir}/libbench/libbench.a @FLIBS@
