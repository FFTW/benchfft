
PRG=doit

INCLUDES = $(INCLBENCH)
PKG=FFTReal
TARFILE=$(PKG).tar.gz
EXTRA_DIST = $(TARFILE)
LOCAL_COPIES = FFTReal.cc FFTReal.h

doit_SOURCES=doit.cc FFTReal.cc FFTReal.h
BUILT_SOURCES = FFTReal.cc FFTReal.h
doit_LDADD=$(LIBBENCH)

include ../Makefile.common

FFTReal.o: FFTReal.h

FFTReal.cc: extract
	@CP@ -f $(PKG)/FFTReal.cpp $@

if DOUBLE
# also replace float with double, as documented in the FFTReal.h header file
FFTReal.h: extract
	sed -e 's/def float.*flt_t/def double flt_t/' $(PKG)/FFTReal.h >$@
else
FFTReal.h: extract
	@CP@ -f $(PKG)/FFTReal.h $@
endif
