if DOUBLE
PRG=doit
endif

PKG = spiral_fft
TARFILE = $(PKG).tar.gz
EXTRA_DIST = $(TARFILE)

INCLUDES = $(INCLBENCH) -I$(PKG)/include

SPIRAL = $(PKG)/lib/libfft.a

doit_SOURCES=doit.c
doit_LDADD=$(SPIRAL) $(LIBBENCH) @FLIBS@
doit.o: $(SPIRAL)

$(SPIRAL): extract
	cd $(PKG) && touch fft_pclprof && ./configure --prefix=`pwd` && $(MAKE) && $(MAKE) install

include ../Makefile.common
