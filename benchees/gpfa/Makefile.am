if HAVE_F77
PRG=doit
noinst_LIBRARIES=libgpfa.a
endif

INCLUDES = $(INCLBENCH)

# these files are clearly apocryphal and we don't use them
EXTRA_DIST=fft3d.f gpfatest.f timer.f gpfft3.f

doit_SOURCES=doit.c
doit_LDADD=libgpfa.a $(LIBBENCH) @FLIBS@

libgpfa_a_SOURCES=gpfa.f gpfa2f.f gpfa3f.f gpfa5f.f setgpfa.f

include ../Makefile.common