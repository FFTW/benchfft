if HAVE_F77
PRG=doit
endif

INCLUDES = $(INCLBENCH)

# these files are clearly apocryphal and we don't use them
EXTRA_DIST=fft3d.f gpfatest.f timer.f gpfft3.f

doit_SOURCES=doit.c gpfa.f gpfa2f.f gpfa3f.f gpfa5f.f setgpfa.f
doit_LDADD=$(LIBBENCH) @FLIBS@

include ../Makefile.common