PRG=doit
PKG=kissfft-131
PKGSUB = $(PKG)/tools
TARFILE=$(PKG).tar.gz
AM_CPPFLAGS = $(INCLBENCH) -I$(PKG) -I$(PKGSUB)
EXTRA_DIST = $(TARFILE)

if SINGLE
AM_CFLAGS=-Dkiss_fft_scalar=float
else
AM_CFLAGS=-Dkiss_fft_scalar=double
endif


doit_SOURCES=doit.c
nodist_doit_SOURCES=$(PKG)/kiss_fft.c \
	$(PKGSUB)/kiss_fftr.c $(PKGSUB)/kiss_fftnd.c $(PKGSUB)/kiss_fftndr.c
doit_LDADD=$(LIBBENCH)

include ../Makefile.common

$(PKG)/kiss_fft.c $(PKG)/kiss_fft.h $(PKGSUB)/kiss_fftr.c $(PKGSUB)/kiss_fftnd.c $(PKGSUB)/kiss_fftndr.c: extract

doit.o: $(PKG)/kiss_fft.h

clean-local::
	rm -rf benchmark extract $(PKG) *.o $(PRG)
