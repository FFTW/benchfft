if HAVE_SGIMATH
noinst_PROGRAMS=doit
else
endif

INCLUDES = -I${top_srcdir}/libbench

doit_SOURCES=doit.c
doit_LDADD=${top_builddir}/libbench/libbench.a @FLIBS@
