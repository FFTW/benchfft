if HAVE_F77
PRG=doit
noinst_LIBRARIES=libharm.a
endif

AM_CPPFLAGS = $(INCLBENCH)

doit_SOURCES=doit.c
doit_LDADD=libharm.a $(LIBBENCH) @FLIBS@

libharm_a_SOURCES = harm.f harmd.f

include ../Makefile.common
