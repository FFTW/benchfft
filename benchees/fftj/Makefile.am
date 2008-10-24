PKG=fftj-0.2
TARFILE=$(PKG).tar.gz
INCLUDES = $(INCLBENCH) -I$(PKG)
EXTRA_DIST = $(TARFILE)

# choose which version of the library to build
if HAVE_F77
if DOUBLE
if HAVE_SSE2
PRG=doit_fort doit_sse
else
PRG=doit_fort
endif
endif
endif

doit_fort_SOURCES=doit.c
doit_fort_LDADD=libfftj_fort.a $(LIBBENCH) @FLIBS@
doit_fort_CFLAGS=-DFFTJ_FORTRAN

doit_sse_SOURCES=doit.c
doit_sse_LDADD=libfftj_sse.a $(LIBBENCH) @FLIBS@
doit_sse_CFLAGS=-DFFTJ_SSE

# -O3 is specified in the FFTJ Makefile
libfftj_fort.a: extract
	(cd $(PKG) && ($(MAKE) fort FC="$(F77) -O3" LIB="../$@"))

libfftj_sse.a: extract
	(cd $(PKG) && ($(MAKE) sse64 FC="$(F77) -O3" LIB="../$@"))

include ../Makefile.common

clean-local::
	rm -rf benchmark extract $(PKG) *.o $(PRG)
