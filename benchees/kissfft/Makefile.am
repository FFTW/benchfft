PRG=doit
PKG=kiss_fft_04
TARFILE=$(PKG).tar.gz
INCLUDES = $(INCLBENCH) -I$(PKG)
EXTRA_DIST = $(TARFILE)

if SINGLE
AM_CFLAGS=-Dkiss_fft_scalar=float
else
AM_CFLAGS=-Dkiss_fft_scalar=double
endif

doit_SOURCES=doit.c
nodist_doit_SOURCES=$(PKG)/kiss_fft.c
doit_LDADD=$(LIBBENCH)

include ../Makefile.common

$(PKG)/kiss_fft.c $(PKG)/kiss_fft.h: extract-mkdir

doit.o: $(PKG)/kiss_fft.h

clean-local::
	rm -rf benchmark extract-mkdir $(PKG) *.o $(PRG)
